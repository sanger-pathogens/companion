#!/usr/bin/env nextflow

statuslog = Channel.create()

genome_file = file(params.inseq)
ref_annot = file(params.ref_annot)
ref_seq = file(params.ref_seq)
ncrna_models = file(params.ncrna_models)

process contiguate_pseudochromosomes {
    input:
    file genome_file
    file ref_seq

    output:
    file 'pseudo.pseudochr.fasta' into pseudochr_seq
    file 'pseudo.pseudochr.agp' into pseudochr_agp
    file 'pseudo.scafs.fasta' into scaffolds_seq
    file 'pseudo.scafs.agp' into scaffolds_agp
    file 'pseudo.contigs.fasta' into contigs_seq

    """
    PATH=${params.ABACAS_DIR}/:$PATH ${params.ABACAS_DIR}/abacas2.nonparallel.sh \
      ${ref_seq} ${genome_file}
    abacas_combine.lua . pseudo "${params.ABACAS_CHR_PATTERN}" \
      "${params.ABACAS_CHR_PREFIX}" "${params.ABACAS_BIN_CHR}" \
      "${params.ABACAS_SEQ_PREFIX}"
    """
}

// make multiple copies of the sequences (better way to do this?!! XXX)
pseudochr_seq_tRNA = Channel.create()
pseudochr_seq_ncRNA = Channel.create()
pseudochr_seq_exonerate = Channel.create()
pseudochr_seq_augustus = Channel.create()
pseudochr_seq_snap = Channel.create()
pseudochr_seq_augustus_ctg = Channel.create()
pseudochr_seq_make_gaps = Channel.create()
pseudochr_seq_dist = Channel.create()
pseudochr_seq_tmhmm = Channel.create()
pseudochr_seq_orthomcl = Channel.create()
pseudochr_seq_splitsplice = Channel.create()
pseudochr_seq.separate(pseudochr_seq_tRNA, pseudochr_seq_ncRNA,
                       pseudochr_seq_augustus, pseudochr_seq_augustus_ctg,
                       pseudochr_seq_snap, pseudochr_seq_dist,
                       pseudochr_seq_make_gaps, pseudochr_seq_splitsplice,
                       pseudochr_seq_tmhmm, pseudochr_seq_orthomcl,
                       pseudochr_seq_exonerate) { a -> [a, a, a, a, a, a, a, a, a, a, a]}

scaffolds_seq_augustus = Channel.create()
scaffolds_seq_make_gaps = Channel.create()
scaffolds_seq.separate(scaffolds_seq_augustus,
                       scaffolds_seq_make_gaps) { a -> [a, a]}

scaffolds_agp_augustus = Channel.create()
scaffolds_agp_make_gaps = Channel.create()
scaffolds_agp.separate(scaffolds_agp_augustus,
                       scaffolds_agp_make_gaps) { a -> [a, a]}

pseudochr_agp_augustus = Channel.create()
pseudochr_agp_make_gaps = Channel.create()
pseudochr_agp.separate(pseudochr_agp_augustus,
                       pseudochr_agp_make_gaps) { a -> [a, a]}

// TRNA PREDICTION
// ===============

process predict_tRNA {
    input:
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_tRNA

    output:
    file 'aragorn.gff3' into trnas
    stdout into statuslog

    """
    ${params.ARAGORN_DIR}/aragorn -t pseudo.pseudochr.fasta | grep -E -C2 '(nucleotides|Sequence)' > 1
    aragorn_to_gff3.lua < 1 > 2
    ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids 2 > aragorn.gff3
    echo "tRNA finished"
    """
}

// NCRNA PREDICTION
// ================

process cmpress {
    input:
    file ncrna_models

    output:
    file ncrna_models into pressed_model

    """
    cmpress -F ${ncrna_models}
    """
}

ncrna_genome_chunk = pseudochr_seq_ncRNA.splitFasta( by: 3)
process predict_ncRNA {
    input:
    file 'chunk' from ncrna_genome_chunk
    file 'models.cm' from pressed_model.first()

    output:
    file 'cm_out' into cmtblouts
    stdout into statuslog

    """
    cmsearch --tblout cm_out --cut_ga models.cm chunk > /dev/null
    echo "ncRNA finished"
    """
}

process merge_ncrnas {
    input:
    file 'cmtblout' from cmtblouts.collectFile()

    output:
    file 'ncrna.gff3' into ncrnafile
    stdout into statuslog

    """
    infernal_to_gff3.lua < ${cmtblout} > 1
    ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids 1 > ncrna.gff3
    echo "ncRNA merged"
    """
}

// better way to do conditional processes? -> see Gitter chat 2015-01-29
if (params.run_exonerate) {
    process make_ref_peps {
        input:
        file ref_annot
        file ref_seq

        output:
        file 'ref.pep' into ref_pep

        """
        ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids ${ref_annot} > 1
        ${params.GT_DIR}/bin/gt extractfeat -matchdescstart -seqfile ${ref_seq} \
          -join -type CDS -translate -retainids 1 > ref.pep
        """
    }
    exn_prot_chunk = ref_pep.splitFasta( by: 20)
    exn_genome_chunk = pseudochr_seq_exonerate.splitFasta( by: 3)
    process run_exonerate {
        input:
        set file('genome.fasta'), file('prot.fasta') from exn_genome_chunk.spread(exn_prot_chunk)

        output:
        file 'exn_out' into exn_results
        stdout into statuslog

        """
        exonerate -E false --model p2g --showvulgar no --showalignment no \
          --showquerygff no --showtargetgff yes --percent 80 \
          --ryo \"AveragePercentIdentity: %pi\n\" prot.fasta \
           genome.fasta > exn_out
        echo "exonerate finished"
        """
    }
    process make_hints {
        input:
        file 'exnout' from exn_results.collectFile()

        output:
        file 'augustus.hints' into exn_hints
        stdout into statuslog

        """
        ${params.AUGUSTUS_SCRIPTDIR}/exonerate2hints.pl \
          --source=P --maxintronlen=${params.AUGUSTUS_HINTS_MAXINTRONLEN} \
          --in=${exnout} \
          --out=augustus.hints
          echo "hints created"
        """
    }
}

process run_augustus_pseudo {
    input:
    //file 'augustus.hints' from exn_hints
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_augustus

    def hintsfile = ""
    if (params.run_exonerate) {
        hintsfile = "--hintsfile=augustus.hints"
    }

    output:
    file 'augustus.gff3' into augustus_pseudo_gff3
    stdout into statuslog

    """
    export AUGUSTUS_CONFIG_PATH=${params.AUGUSTUS_CONFIG_PATH}
    ${params.AUGUSTUS_DIR}/augustus \
        --species=${params.AUGUSTUS_SPECIES} \
        --stopCodonExcludedFromCDS=false \
        --protein=off --codingseq=off --strand=both \
        --genemodel=${params.AUGUSTUS_GENEMODEL} --gff3=on \
        ${hintsfile} \
        --noInFrameStop=true \
        --extrinsicCfgFile=${params.AUGUSTUS_EXTRINSIC_CFG} \
        pseudo.pseudochr.fasta > augustus.full.tmp
    augustus_to_gff3.lua < augustus.full.tmp \
        | ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids \
        > augustus.full.tmp.2
    augustus_mark_partial.lua augustus.full.tmp.2 > augustus.gff3
    echo "AUGUSTUS finished"
    """
}

process run_augustus_contigs {
    input:
    file 'pseudo.contigs.fasta' from contigs_seq
    file 'pseudo.scaffolds.agp' from scaffolds_agp_augustus
    file 'pseudo.scaffolds.fasta' from scaffolds_seq_augustus
    file 'pseudo.pseudochr.agp' from pseudochr_agp_augustus
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_augustus_ctg

    output:
    file 'augustus.scaf.pseudo.mapped.gff3' into augustus_ctg_gff3
    stdout into statuslog

    """
    export AUGUSTUS_CONFIG_PATH=${params.AUGUSTUS_CONFIG_PATH}
    ${params.AUGUSTUS_DIR}/augustus --species=${params.AUGUSTUS_SPECIES} \
        --stopCodonExcludedFromCDS=false \
        --protein=off --codingseq=off --strand=both --genemodel=partial \
        --gff3=on \
        --noInFrameStop=true \
        pseudo.contigs.fasta > augustus.ctg.tmp
    augustus_to_gff3.lua < augustus.ctg.tmp \
        | ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids \
        > augustus.ctg.tmp.2
    augustus_mark_partial.lua augustus.ctg.tmp.2 > augustus.ctg.gff3

    transform_gff_with_agp.lua \
        augustus.ctg.gff3 \
        pseudo.scaffolds.agp \
        pseudo.contigs.fasta \
        pseudo.scaffolds.fasta | \
        ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids > \
        augustus.ctg.scaf.mapped.gff3
    transform_gff_with_agp.lua \
        augustus.ctg.scaf.mapped.gff3 \
        pseudo.pseudochr.agp \
        pseudo.scaffolds.fasta \
        pseudo.pseudochr.fasta | \
        ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids > \
        augustus.scaf.pseudo.mapped.tmp.gff3
    clean_accessions.lua \
        augustus.scaf.pseudo.mapped.tmp.gff3 | \
        ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids > \
        augustus.scaf.pseudo.mapped.gff3
    """
}

process run_snap {
    input:
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_snap

    output:
    file 'snap.gff3' into snap_gff3
    stdout into statuslog

    """
    ${params.SNAP_DIR}/snap -gff -quiet  ${params.SNAP_MODEL} \
        pseudo.pseudochr.fasta > snap.tmp
    snap_gff_to_gff3.lua snap.tmp > snap.gff3
    """
}

process integrate_genemodels {
    input:
    file 'augustus.full.gff3' from augustus_pseudo_gff3
    file 'augustus.ctg.gff3' from augustus_ctg_gff3
    file 'snap.full.gff3' from snap_gff3

    output:
    file 'integrated.fixed.sorted.gff3' into integrated_gff3

    """
    unset GT_RETAINIDS
    ${params.GT_DIR}/bin/gt gff3 -fixregionboundaries -retainids no -sort -tidy \
        augustus.full.gff3 augustus.ctg.gff3 snap.full.gff3 \
        > merged.gff3
    export GT_RETAINIDS=yes

    # choose best gene model for overlapping region
    integrate_gene_calls.lua merged.gff3 | \
        ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids \
        > integrated.gff3

    # make this parameterisable
    fix_polycistrons.lua integrated.gff3 > integrated.fixed.gff3

    # make sure final output is sorted
    ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids \
      integrated.fixed.gff3 > integrated.fixed.sorted.gff3
    """
}

process merge_structural {
    input:
    file 'ncrna.gff3' from ncrnafile
    file 'trna.gff3' from trnas
    file 'integrated.gff3' from integrated_gff3

    output:
    file 'structural.full.gff3' into genemodels_gff3

    """
    ${params.GT_DIR}/bin/gt gff3 -sort -tidy ncrna.gff3 trna.gff3 integrated.gff3 \
        > structural.full.gff3
    """
}

process add_gap_features {
    input:
    file 'merged_in.gff3' from genemodels_gff3
    file 'pseudo.scaffolds.agp' from scaffolds_agp_make_gaps
    file 'pseudo.scaffolds.fasta' from scaffolds_seq_make_gaps
    file 'pseudo.pseudochr.agp' from pseudochr_agp_make_gaps
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_make_gaps

    output:
    file 'merged_out.gff3' into genemodels_with_gaps_gff3

    """
    make_contig_features_from_agp.lua pseudo.scaffolds.agp "within scaffold" | \
      ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids > contigs.gff3

    transform_gff_with_agp.lua contigs.gff3 \
      pseudo.pseudochr.agp pseudo.scaffolds.fasta pseudo.pseudochr.fasta "between scaffolds" | \
      ${params.GT_DIR}/bin/gt gff3 -sort -tidy -retainids > contigs2.gff3

    ${params.GT_DIR}/bin/gt merge -force -o merged_out.gff3 \
      merged_in.gff3 contigs2.gff3
    """
}

process split_splice_models_at_gaps {
    input:
    file 'input.gff3' from genemodels_with_gaps_gff3
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_splitsplice

    output:
    file 'merged_out.gff3' into genemodels_with_gaps_split_gff3

    """
    # sort
    ${params.GT_DIR}/bin/gt gff3 -sort -retainids -tidy < input.gff3 > tmp2

    # split genes at inter-scaffold gaps
    split_genes_at_gaps.lua tmp2 | \
      ${params.GT_DIR}/bin/gt gff3 -sort -retainids -tidy > tmp3

    # splice genes at inter-contig gaps, get rid of short partials
    splice_genes_at_gaps.lua tmp3 | \
      ${params.GT_DIR}/bin/gt gff3 -sort -retainids -tidy | \
      ${params.GT_DIR}/bin/gt select \
          -rule_files ${params.FILTER_SHORT_PARTIALS_RULE} \
      > tmp4

    # get rid of genes still having stop codons
    filter_genes_with_stop_codons.lua \
      tmp4 pseudo.pseudochr.fasta | \
      ${params.GT_DIR}/bin/gt gff3 -sort -retainids -tidy > merged_out.gff3
    """
}

process add_polypeptides {
    input:
    file 'input.gff3' from genemodels_with_gaps_split_gff3

    output:
    file 'output.gff3' into genemodels_with_polypeptides_gff3

    """
    create_polypeptides.lua input.gff3 ${params.GENOME_PREFIX} \
        "${params.CHR_PATTERN}" > output.gff3
    """
}

process run_tmhmm {
    input:
    file 'input.gff3' from genemodels_with_polypeptides_gff3
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_tmhmm

    output:
    file 'output.gff3' into genemodels_with_tmhmm_gff3

    """
    ${params.GT_DIR}/bin/gt extractfeat -type CDS -join -translate \
        -seqfile pseudo.pseudochr.fasta -matchdescstart \
        input.gff3 > proteins.fas
    ${params.TMHMM_DIR}/tmhmm --noplot < proteins.fas > tmhmm.out
    tmhmm_to_gff3.lua tmhmm.out input.gff3 | \
        ${params.GT_DIR}/bin/gt gff3 -sort -retainids -tidy > output.gff3
    """
}

process get_proteins_for_orthomcl {
    input:
    file 'input.gff3' from genemodels_with_tmhmm_gff3
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_orthomcl

    output:
    file 'input.gff3' into genemodels_orthomcl_gff3
    file 'proteins.fas' into proteins_orthomcl_gff3

    """
    ${params.GT_DIR}/bin/gt extractfeat -type CDS -translate -join -retainids \
        -seqfile pseudo.pseudochr.fasta -matchdescstart < input.gff3 \
        | truncate_header.lua > proteins.fas


    ${params.GT_DIR}/bin/gt extractfeat -type pseudogenic_exon -translate \
        -join -retainids \
        -seqfile pseudo.pseudochr.fasta -matchdescstart < input.gff3 \
        | truncate_header.lua >> proteins.fas
    """
}

pepfiles = Channel.from(params.OMCL_PEPFILES)
process make_ref_input_for_orthomcl {
    tag { shortname }

    input:
    set shortname, pepfile from pepfiles

    output:
    file 'out.gg' into gg_file
    file 'out.map' into mapfile
    file 'shortname' into shortname
    file 'mapped.fasta' into mapped_fasta

    """
    truncate_header.lua < ${pepfile} > pepfile.trunc
    map_protein_names.lua ${shortname} pepfile.trunc out.map > mapped.fasta
    make_gg_line.lua ${shortname} mapped.fasta > out.gg
    echo "${shortname}" > shortname
    """
}

process make_target_input_for_orthomcl {
    input:
    file 'pepfile.fas' from proteins_orthomcl_gff3

    output:
    file 'out.gg' into gg_file_ref
    file 'out.map' into mapfile_ref
    file 'shortname' into shortname_ref
    file 'mapped.fasta' into mapped_fasta_ref

    """
    truncate_header.lua < pepfile.fas > pepfile.trunc
    map_protein_names.lua ${params.GENOME_PREFIX} pepfile.trunc out.map > mapped.fasta
    make_gg_line.lua ${params.GENOME_PREFIX} mapped.fasta > out.gg
    echo "${params.GENOME_PREFIX}" > shortname
    """
}
full_gg = gg_file.mix(gg_file_ref).collectFile()
full_shortnames = shortname.mix(shortname_ref).collectFile()
full_mapped_fasta = mapped_fasta.mix(mapped_fasta_ref).collectFile()
full_mapfile = mapfile.mix(mapfile_ref).collectFile()

process blast_for_orthomcl {
    cpus 8
    echo true

    input:
    file 'mapped.fasta' from full_mapped_fasta

    output:
    file 'blastout' into orthomcl_blastout

    """
    formatdb -i mapped.fasta
    blastall -p blastp -F 'm S' -a 10 -e 1e-5 -d mapped.fasta \
      -m 8 -i mapped.fasta > blastout
    """
}

process run_orthomcl {
    echo true

    input:
    file 'blastout' from orthomcl_blastout
    file 'ggfile' from full_gg

    output:
    file 'orthomcl_out' into orthomcl_cluster_out

    """
    ${params.OMCL_DIR}/orthomcl.pl --mode 3 \
      --blast_file blastout \
      --gg_file ggfile

    ORTHOMCL_OUTFILE =`find . -mindepth 1 -name all_orthomcl.out`
    if [ -e $ORTHOMCL_OUTFILE ]; then
      cp $ORTHOMCL_OUTFILE orthomcl_out
    fi
    """
}

// ========  REPORTING ========

process report {
    input:
    file 'models.gff3' from genemodels_orthomcl_gff3

    output:
    stdout into result

    """
    ls -HAl models.gff3
    ${params.GT_DIR}/bin/gt stat models.gff3
 #   cat models.gff3
    """
}

result.subscribe {
    println it
}

statuslog.subscribe {
   // println "log:  " + it
}