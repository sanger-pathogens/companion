#!/usr/bin/env nextflow

statuslog = Channel.create()

genome_file = file(params.inseq)
ref_annot = file(params.ref_annot)
ref_seq = file(params.ref_seq)
go_obo = file(params.GO_OBO)
ncrna_models = file(params.NCRNA_MODELS)

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
    abacas2.nonparallel.sh \
      ${ref_seq} ${genome_file}
    abacas_combine.lua . pseudo "${params.ABACAS_CHR_PATTERN}" \
      "${params.ABACAS_CHR_PREFIX}" "${params.ABACAS_BIN_CHR}" \
      "${params.ABACAS_SEQ_PREFIX}"
    """
}

// fork important output streams
pseudochr_seq_tRNA = Channel.create()
pseudochr_seq_ncRNA = Channel.create()
pseudochr_seq_ratt = Channel.create()
pseudochr_seq_exonerate = Channel.create()
pseudochr_seq_augustus = Channel.create()
pseudochr_seq_snap = Channel.create()
pseudochr_seq_augustus_ctg = Channel.create()
pseudochr_seq_make_gaps = Channel.create()
pseudochr_seq_dist = Channel.create()
pseudochr_seq_tmhmm = Channel.create()
pseudochr_seq_orthomcl = Channel.create()
pseudochr_seq_splitsplice = Channel.create()
pseudochr_seq.into(pseudochr_seq_tRNA, pseudochr_seq_ncRNA, pseudochr_seq_ratt,
                   pseudochr_seq_augustus, pseudochr_seq_augustus_ctg,
                   pseudochr_seq_snap, pseudochr_seq_dist,
                   pseudochr_seq_make_gaps, pseudochr_seq_splitsplice,
                   pseudochr_seq_tmhmm, pseudochr_seq_orthomcl,
                   pseudochr_seq_exonerate)

scaffolds_seq_augustus = Channel.create()
scaffolds_seq_make_gaps = Channel.create()
scaffolds_seq.into(scaffolds_seq_augustus, scaffolds_seq_make_gaps)

scaffolds_agp_augustus = Channel.create()
scaffolds_agp_make_gaps = Channel.create()
scaffolds_agp.into(scaffolds_agp_augustus, scaffolds_agp_make_gaps)

pseudochr_agp_augustus = Channel.create()
pseudochr_agp_make_gaps = Channel.create()
pseudochr_agp.into(pseudochr_agp_augustus, pseudochr_agp_make_gaps)

// TRNA PREDICTION
// ===============

process predict_tRNA {
    input:
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_tRNA

    output:
    file 'aragorn.gff3' into trnas
    stdout into statuslog

    """
    aragorn -t pseudo.pseudochr.fasta \
        | grep -E -C2 '(nucleotides|Sequence)' > 1
    aragorn_to_gff3.lua < 1 > 2
    gt gff3 -sort -tidy -retainids 2 > aragorn.gff3
    echo "tRNA finished"
    """
}

// NCRNA PREDICTION
// ================

ncrna_genome_chunk = pseudochr_seq_ncRNA.splitFasta( by: 3)
process predict_ncRNA {
    input:
    file 'chunk' from ncrna_genome_chunk
    val ncrna_models

    output:
    file 'cm_out' into cmtblouts
    stdout into statuslog

    """
    cmpress -F ${ncrna_models}
    cmsearch --tblout cm_out --cut_ga ${ncrna_models} chunk > /dev/null
    echo "ncRNA finished"
    """
}

process merge_ncrnas {
    cache 'deep'

    input:
    file 'cmtblout' from cmtblouts.collectFile()

    output:
    file 'ncrna.gff3' into ncrnafile
    stdout into statuslog

    """
    infernal_to_gff3.lua < ${cmtblout} > 1
    gt gff3 -sort -tidy -retainids 1 > ncrna.gff3
    echo "ncRNA merged"
    """
}

// PROTEIN-DNA ALIGNMENT
// =====================
if (params.run_exonerate) {
    process make_ref_peps {
        cache 'deep'

        input:
        file ref_annot
        file ref_seq

        output:
        file 'ref.pep' into ref_pep

        """
        gt gff3 -sort -tidy -retainids ${ref_annot} > 1
        gt extractfeat -matchdescstart -seqfile ${ref_seq} \
          -join -type CDS -translate -retainids 1 > ref.pep
        """
    }
    exn_prot_chunk = ref_pep.splitFasta( by: 20)
    exn_genome_chunk = pseudochr_seq_exonerate.splitFasta( by: 3)
    process run_exonerate {
        cache 'deep'

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
        cache 'deep'

        input:
        file 'exnout' from exn_results.collectFile()

        output:
        set val("--hintsfile=augustus.hints"), file('augustus.hints') into exn_hints
        stdout into statuslog

        """
        ${params.AUGUSTUS_SCRIPTDIR}/exonerate2hints.pl \
          --source=P --maxintronlen=${params.AUGUSTUS_HINTS_MAXINTRONLEN} \
          --in=${exnout} \
          --out=augustus.hints
          echo "hints created"
        """
    }
} else {
    exn_hints = Channel.just(["",file("/dev/null")])
}

// RATT
// ====

process ratt_make_ref_embl {
    input:
    file ref_annot
    file ref_seq
    val go_obo

    output:
    file '*.embl' into ref_embl
    stdout into statuslog

    """
    # split away
    gt inlineseq_split -seqfile /dev/null -gff3file ref_without_seq.gff3 ${ref_annot}
    gt inlineseq_add -seqfile ${ref_seq} -matchdescstart ref_without_seq.gff3 > ref_with_seq.gff3
    gff3_to_embl.lua ref_with_seq.gff3 ${go_obo} Foo
    """
}

process run_ratt {
    input:
    file 'in*.embl' from ref_embl
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_ratt

    output:
    file 'Out*.final.embl' into ratt_result
    file 'Out*.Report.txt' into ratt_reports
    stdout into statuslog

    """
    start.ratt.sh . pseudo.pseudochr.fasta Out ${params.RATT_TRANSFER_TYPE}
    """
}

process ratt_to_gff3 {
    input:
    file 'in*.embl' from ratt_result
    file 'in*.report' from ratt_reports

    output:
    file 'ratt.gff3' into ratt_gff3
    stdout into statuslog

    """
    ratt_embl_to_gff3.lua in*.embl | \
      gt gff3 -sort -retainids -tidy > \
      ratt.tmp.gff3
    ratt_remove_problematic.lua ratt.tmp.gff3 in*report | \
      gt gff3 -sort -retainids -tidy > \
      ratt.gff3
    """
}

// AUGUSTUS
// ================

process run_augustus_pseudo {
    cache 'deep'

    input:
    set val(hintsline), file('augustus.hints') from exn_hints
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_augustus

    output:
    file 'augustus.gff3' into augustus_pseudo_gff3
    stdout into statuslog

    """
    augustus \
        --species=${params.AUGUSTUS_SPECIES} \
        --stopCodonExcludedFromCDS=false \
        --protein=off --codingseq=off --strand=both \
        --genemodel=${params.AUGUSTUS_GENEMODEL} --gff3=on \
        ${hintsline} \
        --noInFrameStop=true \
        --extrinsicCfgFile=${params.AUGUSTUS_EXTRINSIC_CFG} \
        pseudo.pseudochr.fasta > augustus.full.tmp
    augustus_to_gff3.lua < augustus.full.tmp \
        | gt gff3 -sort -tidy -retainids \
        > augustus.full.tmp.2
    augustus_mark_partial.lua augustus.full.tmp.2 > augustus.gff3
    echo "AUGUSTUS finished"
    """
}

process run_augustus_contigs {
    cache 'deep'

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
    augustus --species=${params.AUGUSTUS_SPECIES} \
        --stopCodonExcludedFromCDS=false \
        --protein=off --codingseq=off --strand=both --genemodel=partial \
        --gff3=on \
        --noInFrameStop=true \
        pseudo.contigs.fasta > augustus.ctg.tmp
    augustus_to_gff3.lua < augustus.ctg.tmp \
        | gt gff3 -sort -tidy -retainids \
        > augustus.ctg.tmp.2
    augustus_mark_partial.lua augustus.ctg.tmp.2 > augustus.ctg.gff3

    transform_gff_with_agp.lua \
        augustus.ctg.gff3 \
        pseudo.scaffolds.agp \
        pseudo.contigs.fasta \
        pseudo.scaffolds.fasta | \
        gt gff3 -sort -tidy -retainids > \
        augustus.ctg.scaf.mapped.gff3
    transform_gff_with_agp.lua \
        augustus.ctg.scaf.mapped.gff3 \
        pseudo.pseudochr.agp \
        pseudo.scaffolds.fasta \
        pseudo.pseudochr.fasta | \
        gt gff3 -sort -tidy -retainids > \
        augustus.scaf.pseudo.mapped.tmp.gff3
    clean_accessions.lua \
        augustus.scaf.pseudo.mapped.tmp.gff3 | \
        gt gff3 -sort -tidy -retainids > \
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
    snap -gff -quiet  ${SNAP_MODEL} \
        pseudo.pseudochr.fasta > snap.tmp
    snap_gff_to_gff3.lua snap.tmp > snap.gff3
    """
}

process integrate_genemodels {
    cache 'deep'

    input:
    file 'augustus.full.gff3' from augustus_pseudo_gff3
    file 'augustus.ctg.gff3' from augustus_ctg_gff3
    file 'snap.full.gff3' from snap_gff3
    file 'ratt.full.gff3' from ratt_gff3

    output:
    file 'integrated.fixed.sorted.gff3' into integrated_gff3

    """
    unset GT_RETAINIDS
    gt gff3 -fixregionboundaries -retainids no -sort -tidy \
        augustus.full.gff3 augustus.ctg.gff3 snap.full.gff3 ratt.full.gff3 \
        > merged.gff3
    export GT_RETAINIDS=yes

    # choose best gene model for overlapping region
    integrate_gene_calls.lua merged.gff3 | \
        gt gff3 -sort -tidy -retainids \
        > integrated.gff3

    # TODO: make this parameterisable
    fix_polycistrons.lua integrated.gff3 > integrated.fixed.gff3

    # make sure final output is sorted
    gt gff3 -sort -tidy -retainids \
      integrated.fixed.gff3 > integrated.fixed.sorted.gff3
    """
}

// MERGE ALL GENES TO FINAL SET AND CLEANUP
// ========================================

process merge_structural {
    cache 'deep'

    input:
    file 'ncrna.gff3' from ncrnafile
 //   file 'trna.gff3' from trnas
    file 'integrated.gff3' from integrated_gff3

    output:
    file 'structural.full.gff3' into genemodels_gff3

    // removed trna.gff3
    """
    gt gff3 -sort -tidy ncrna.gff3 integrated.gff3 \
        > structural.full.gff3
    """
}

// ADD GAPS AND ADJUST GENES NOT TO SPAN GAPS
// ==========================================

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
      gt gff3 -sort -tidy -retainids > contigs.gff3

    transform_gff_with_agp.lua contigs.gff3 \
      pseudo.pseudochr.agp pseudo.scaffolds.fasta \
      pseudo.pseudochr.fasta "between scaffolds" | \
      gt gff3 -sort -tidy -retainids > contigs2.gff3

    gt merge -force -o merged_out.gff3 \
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
    gt gff3 -sort -retainids -tidy < input.gff3 > tmp2

    # split genes at inter-scaffold gaps
    split_genes_at_gaps.lua tmp2 | \
      gt gff3 -sort -retainids -tidy > tmp3

    # splice genes at inter-contig gaps, get rid of short partials
    splice_genes_at_gaps.lua tmp3 | \
      gt gff3 -sort -retainids -tidy | \
      gt select \
          -rule_files ${FILTER_SHORT_PARTIALS_RULE} \
      > tmp4

    # get rid of genes still having stop codons
    filter_genes_with_stop_codons.lua \
      tmp4 pseudo.pseudochr.fasta | \
      gt gff3 -sort -retainids -tidy > merged_out.gff3
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

// TMHMM
// =====

// disabled for now for license reasons
/* process run_tmhmm {
    input:
    file 'input.gff3' from genemodels_with_polypeptides_gff3
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_tmhmm

    output:
    file 'output.gff3' into genemodels_with_tmhmm_gff3

    """
    gt extractfeat -type CDS -join -translate \
        -seqfile pseudo.pseudochr.fasta -matchdescstart \
        input.gff3 > proteins.fas
    tmhmm --noplot < proteins.fas > tmhmm.out
    tmhmm_to_gff3.lua tmhmm.out input.gff3 | \
        gt gff3 -sort -retainids -tidy > output.gff3
    """
} */

// ORTHOMCL AND FUNCTIONAL TRANSFER
// ================================

genemodels_for_omcl_proteins = Channel.create()
genemodels_for_omcl_annot = Channel.create()
genemodels_with_polypeptides_gff3.into(genemodels_for_omcl_proteins,
                                       genemodels_for_omcl_annot)

process get_proteins_for_orthomcl {
    input:
    file 'input.gff3' from genemodels_for_omcl_proteins
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_orthomcl

    output:
    file 'input.gff3' into genemodels_orthomcl_gff3
    file 'proteins.fas' into proteins_target

    """
    gt extractfeat -type CDS -translate -join -retainids \
        -seqfile pseudo.pseudochr.fasta -matchdescstart < input.gff3 \
        | truncate_header.lua > proteins.fas


    gt extractfeat -type pseudogenic_exon -translate \
        -join -retainids \
        -seqfile pseudo.pseudochr.fasta -matchdescstart < input.gff3 \
        | truncate_header.lua >> proteins.fas
    """
}
proteins_orthomcl = Channel.create()
proteins_pfam = Channel.create()

proteins_target.into(proteins_orthomcl, proteins_pfam)

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
    file 'pepfile.fas' from proteins_orthomcl

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

full_mapped_fasta_for_index = Channel.create()
full_mapped_fasta_for_query = Channel.create()
full_mapped_fasta.into(full_mapped_fasta_for_index, full_mapped_fasta_for_query)

process blast_for_orthomcl_formatdb {
    cache 'deep'

    input:
    file 'mapped.fasta' from full_mapped_fasta_for_index

    output:
    file 'mapped.fasta' into full_mapped_fasta_indexed
    file 'mapped.fasta.phr' into full_mapped_fasta_indexed_phr
    file 'mapped.fasta.psq' into full_mapped_fasta_indexed_psq
    file 'mapped.fasta.pin' into full_mapped_fasta_indexed_pin

    """
    formatdb -i mapped.fasta
    """
}

proteins_orthomcl_blast_chunk = full_mapped_fasta_for_query.splitFasta( by: 10)
process blast_for_orthomcl {
    cache 'deep'

    input:
    file 'mapped_chunk.fasta' from proteins_orthomcl_blast_chunk
    file 'mapped.fasta' from full_mapped_fasta_indexed.first()
    file 'mapped.fasta.phr' from full_mapped_fasta_indexed_phr.first()
    file 'mapped.fasta.psq' from full_mapped_fasta_indexed_psq.first()
    file 'mapped.fasta.pin' from full_mapped_fasta_indexed_pin.first()

    output:
    file 'blastout' into orthomcl_blastout

    """
    blastall -p blastp -F 'm S' -e 1e-5 -d mapped.fasta \
      -m 8 -i mapped_chunk.fasta > blastout
    """
}

process run_orthomcl {
    cache 'deep'

    input:
    file 'blastout' from orthomcl_blastout.collectFile()
    file 'ggfile' from full_gg

    output:
    file 'orthomcl_out' into orthomcl_cluster_out

    """
    orthomcl.pl --inflation 1.5 --mode 3 \
      --blast_file blastout \
      --gg_file ggfile
    cp `find . -mindepth 1 -name all_orthomcl.out` orthomcl_out
    """
}

process annotate_orthologs {
    cache 'deep'

    input:
    file 'orthomcl_out' from orthomcl_cluster_out
    file 'mapfile' from full_mapfile
    file 'input.gff3' from genemodels_for_omcl_annot

    output:
    file 'with_func.gff3' into gff3_with_ortho_transferred
    file 'orthomcl.gaf' into gaf_with_ortho_transferred

    """
    # annotate GFF with ortholog clusters and members
    map_clusters_gff.lua input.gff3 orthomcl_out mapfile > with_clusters.gff3

    # transfer functional annotation from orthologs
    transfer_annotations_from_gff.lua with_clusters.gff3 \
        ${params.OMCL_GFFFILE} > with_func.gff3

    # transfer GOs from orthologs GAF
    transfer_annotations_from_gaf.lua with_func.gff3 \
        ${params.OMCL_GAFFILE} ${params.DB_ID} \
        ${params.TAXON_ID} > orthomcl.gaf
    """
}

// PFAM DOMAIN ANNOTATION
// ======================

proteins_pfam_chunk = proteins_pfam.splitFasta( by: 10)
process run_pfam {
    cache 'deep'

    input:
    file 'proteins.fas' from proteins_pfam_chunk

    output:
    file 'pfamout' into pfam_output

    """
    hmmscan --domtblout pfamout --cut_ga --noali --cpu 2 ${PFAM} proteins.fas
    """
}

process pfam_to_gff3 {
    cache 'deep'

    input:
    file 'pfam.out' from pfam_output.collectFile()

    output:
    file 'pfam.gff3' into pfam_gff3

    """
    pfam_to_gff3.lua ${PFAM2GO} < pfam.out | gt gff3 -sort -tidy -retainids > pfam.gff3
    """
}

process annotate_pfam {
    cache 'deep'

    input:
    file 'pfam.gff3' from pfam_gff3
    file 'with_ortho.gff3' from gff3_with_ortho_transferred
    file 'with_ortho.gaf' from gaf_with_ortho_transferred
    val go_obo

    output:
    file 'with_pfam.gff3' into gff3_with_pfam
    file 'with_pfam.gaf' into gaf_with_pfam

    """
    iproscan_gff3_merge.lua with_ortho.gff3 pfam.gff3 | \
      gt gff3 -tidy -sort -retainids > with_pfam.gff3
    iproscan_gaf_merge.lua with_pfam.gff3 pfam.gff3 \
      ${params.DB_ID} ${params.TAXON_ID} ${go_obo} > my.with_pfam.gaf
    cat with_ortho.gaf my.with_pfam.gaf > with_pfam.gaf
    """
}


// REPORTING
// =========

process report {
    input:
    file 'with_pfam.gff3' from gff3_with_pfam
    file 'with_pfam.gaf' from gaf_with_pfam

    output:
    file 'with_pfam.gff3' into result_gff3
    file 'with_pfam.gaf' into result_gaf

    """
    ls -HAl with_pfam.gff3
    """
}

result_gff3.collectFile().subscribe {
    println it
}

result_gaf.collectFile().subscribe {
    println it
}

statuslog.subscribe {
   // println "log:  " + it
}
