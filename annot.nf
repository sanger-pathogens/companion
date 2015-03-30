#!/usr/bin/env nextflow

/*
    Copyright (c) 2014-2015 Sascha Steinbiss <ss34@sanger.ac.uk>
    Copyright (c) 2014-2015 Genome Research Ltd

    Permission to use, copy, modify, and distribute this software for any
    purpose with or without fee is hereby granted, provided that the above
    copyright notice and this permission notice appear in all copies.

    THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
    WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
    ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
    WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
    ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
    OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

// TODO check for required parameters!

genome_file = file(params.inseq)
ref_annot = file(params.ref_dir + "/" + params.ref_species + "/annotation.gff3")
ref_seq = file(params.ref_dir + "/" + params.ref_species + "/genome.fasta")
ref_chr = file(params.ref_dir + "/" + params.ref_species + "/chromosomes.fasta")
go_obo = file(params.GO_OBO)
ncrna_models = file(params.NCRNA_MODELS)
extrinsic_cfg = file(params.AUGUSTUS_EXTRINSIC_CFG)
omcl_gfffile = file(params.ref_dir + "/" + params.ref_species + "/annotation.gff3")
omcl_gaffile = file(params.ref_dir + "/" + params.ref_species + "/go.gaf")
omcl_pepfile = file(params.ref_dir + "/" + params.ref_species + "/proteins.fasta")

// PSEUDOCHROMOSOME CONTIGUATION
// =============================

if (params.do_contiguation) {
    process contiguate_pseudochromosomes {
        input:
        file genome_file
        file ref_chr
        val params.ABACAS_CHR_PATTERN
        val params.ABACAS_CHR_PREFIX
        val params.ABACAS_BIN_CHR
        val params.ABACAS_SEQ_PREFIX

        output:
        file 'pseudo.pseudochr.fasta' into pseudochr_seq
        file 'pseudo.pseudochr.agp' into pseudochr_agp
        file 'pseudo.scafs.fasta' into scaffolds_seq
        file 'pseudo.scafs.agp' into scaffolds_agp
        file 'pseudo.contigs.fasta' into contigs_seq
        file 'ref_target_mapping.txt' into ref_target_mapping

        """
        abacas2.nonparallel.sh \
          ${ref_chr} ${genome_file} 500 85 0 3000
        abacas_combine.lua . pseudo "${params.ABACAS_CHR_PATTERN}" \
          "${params.ABACAS_CHR_PREFIX}" "${params.ABACAS_BIN_CHR}" \
          "${params.ABACAS_SEQ_PREFIX}"
        """
    }
} else {
    process prepare_noncontiguated_input {
        input:
        file genome_file

        output:
        file 'pseudo.pseudochr.fasta' into pseudochr_seq
        file 'pseudo.pseudochr.agp' into pseudochr_agp
        file 'pseudo.scafs.fasta' into scaffolds_seq
        file 'pseudo.scafs.agp' into scaffolds_agp
        file 'pseudo.contigs.fasta' into contigs_seq

        """
        no_abacas_prepare.lua ${genome_file} pseudo
        """
    }
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
pseudochr_seq_make_dist_1 = Channel.create()
pseudochr_seq_make_dist_2 = Channel.create()
pseudochr_seq_tmhmm = Channel.create()
pseudochr_seq_orthomcl = Channel.create()
pseudochr_seq_splitsplice = Channel.create()
pseudochr_seq.into(pseudochr_seq_tRNA, pseudochr_seq_ncRNA, pseudochr_seq_ratt,
                   pseudochr_seq_augustus, pseudochr_seq_augustus_ctg,
                   pseudochr_seq_snap, pseudochr_seq_make_gaps,
                   pseudochr_seq_make_dist_1, pseudochr_seq_make_dist_2,
                   pseudochr_seq_splitsplice,
                   pseudochr_seq_tmhmm, pseudochr_seq_orthomcl,
                   pseudochr_seq_exonerate)

scaffolds_seq_augustus = Channel.create()
scaffolds_seq_make_gaps = Channel.create()
scaffolds_seq_make_dist_1 = Channel.create()
scaffolds_seq_make_dist_2 = Channel.create()
scaffolds_seq.into(scaffolds_seq_augustus, scaffolds_seq_make_gaps,
                   scaffolds_seq_make_dist_1, scaffolds_seq_make_dist_2)

scaffolds_agp_augustus = Channel.create()
scaffolds_agp_make_gaps = Channel.create()
scaffolds_agp_make_dist = Channel.create()
scaffolds_agp.into(scaffolds_agp_augustus, scaffolds_agp_make_gaps,
                   scaffolds_agp_make_dist)

pseudochr_agp_augustus = Channel.create()
pseudochr_agp_make_gaps = Channel.create()
pseudochr_agp_make_dist = Channel.create()
pseudochr_agp.into(pseudochr_agp_augustus, pseudochr_agp_make_gaps,
                   pseudochr_agp_make_dist)

// TRNA PREDICTION
// ===============

process predict_tRNA {
    input:
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_tRNA

    output:
    file 'aragorn.gff3' into trnas

    """
    aragorn -t pseudo.pseudochr.fasta > out
    grep -E -C2 '(nucleotides|Sequence)' out > 1
    aragorn_to_gff3.lua < 1 > 2
    gt gff3 -sort -tidy -retainids 2 > aragorn.gff3
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

    """
    cp ${ncrna_models} .
    cmpress -F *.cm
    cmsearch --tblout cm_out --cut_ga ${ncrna_models} chunk > /dev/null
    """
}

process merge_ncrnas {
    cache 'deep'

    input:
    file 'cmtblout' from cmtblouts.collectFile()

    output:
    file 'ncrna.gff3' into ncrnafile

    """
    infernal_to_gff3.lua < ${cmtblout} > 1
    gt gff3 -sort -tidy -retainids 1 > ncrna.gff3
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

        """
        exonerate -E false --model p2g --showvulgar no --showalignment no \
          --showquerygff no --showtargetgff yes --percent 80 \
          --ryo \"AveragePercentIdentity: %pi\n\" prot.fasta \
           genome.fasta > exn_out
        """
    }
    process make_hints {
        cache 'deep'

        input:
        file 'exnout' from exn_results.collectFile()

        output:
        set val(outline), file('augustus.hints') into exn_hints

        script:
        outline = "--hintsfile=augustus.hints"
        """
        exonerate2hints.pl \
          --source=P --maxintronlen=${params.AUGUSTUS_HINTS_MAXINTRONLEN} \
          --in=${exnout} \
          --out=augustus.hints
        """
    }
} else {
    exn_hints = Channel.just(["",""])
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

    """
    # make sure GFF3 contains sequence
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
    val extrinsic_cfg

    output:
    file 'augustus.gff3' into augustus_pseudo_gff3

    """
    augustus \
        --species=${params.AUGUSTUS_SPECIES} \
        --stopCodonExcludedFromCDS=false \
        --protein=off --codingseq=off --strand=both \
        --genemodel=${params.AUGUSTUS_GENEMODEL} --gff3=on \
        ${hintsline} \
        --noInFrameStop=true \
        --extrinsicCfgFile=${extrinsic_cfg} \
        pseudo.pseudochr.fasta > augustus.full.tmp
    augustus_to_gff3.lua < augustus.full.tmp \
        | gt gff3 -sort -tidy -retainids \
        | gt select -mingenescore ${params.AUGUSTUS_SCORE_THRESHOLD} \
        > augustus.full.tmp.2
    augustus_mark_partial.lua augustus.full.tmp.2 > augustus.gff3
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

    """
    augustus --species=${params.AUGUSTUS_SPECIES} \
        --stopCodonExcludedFromCDS=false \
        --protein=off --codingseq=off --strand=both --genemodel=partial \
        --gff3=on \
        --noInFrameStop=true \
        pseudo.contigs.fasta > augustus.ctg.tmp && \
    augustus_to_gff3.lua < augustus.ctg.tmp \
        | gt gff3 -sort -tidy -retainids \
        | gt select -mingenescore ${params.AUGUSTUS_SCORE_THRESHOLD} \
        > augustus.ctg.tmp.2 && \
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
        augustus.scaf.pseudo.mapped.gff3
    #clean_accessions.lua \
    #    augustus.scaf.pseudo.mapped.tmp.gff3 | \
    #    gt gff3 -sort -tidy -retainids > \
    #    augustus.scaf.pseudo.mapped.gff3
    """
}

process run_snap {
    input:
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_snap

    output:
    file 'snap.gff3' into snap_gff3

    """
    snap -gff -quiet  ${SNAP_MODEL} \
        pseudo.pseudochr.fasta > snap.tmp
    snap_gff_to_gff3.lua snap.tmp > snap.tmp.2
    gt gff3 -sort -tidy -retainids snap.tmp.2 > snap.gff3
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
        > merged.pre.gff3
    export GT_RETAINIDS=yes

    # avoid huge gene clusters
    gt select -maxgenelength ${params.MAX_GENE_LENGTH} merged.pre.gff3 > merged.gff3

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

process remove_exons {
    input:
    file 'integrated.gff3' from integrated_gff3

    output:
    file 'integrated_clean.gff3' into integrated_gff3_clean

    """
    remove_exons.lua integrated.gff3 > integrated_clean.gff3
    """
}

// MERGE ALL GENES TO FINAL SET AND CLEANUP
// ========================================

process merge_structural {
    cache 'deep'

    input:
    file 'ncrna.gff3' from ncrnafile
    file 'trna.gff3' from trnas
    file 'integrated.gff3' from integrated_gff3_clean

    output:
    file 'structural.full.gff3' into genemodels_gff3

    """
    gt gff3 -sort -tidy ncrna.gff3 integrated.gff3 trna.gff3 \
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
    set -ev
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
    val params.CHR_PATTERN

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
refcomp_protein_in = Channel.create()
proteins_output = Channel.create()
proteins_target.into(proteins_orthomcl, proteins_pfam, refcomp_protein_in,
                     proteins_output)

process make_ref_input_for_orthomcl {
    input:
    val omcl_pepfile
    val params.ref_species

    output:
    file 'out.gg' into gg_file
    file 'shortname' into shortname
    file 'mapped.fasta' into mapped_fasta

    script:
    """
    truncate_header.lua < ${omcl_pepfile} > pepfile.trunc
    ln -s pepfile.trunc mapped.fasta
    make_gg_line.lua ${params.ref_species} mapped.fasta > out.gg
    echo "${params.ref_species}" > shortname
    """
}

process make_target_input_for_orthomcl {
    input:
    file 'pepfile.fas' from proteins_orthomcl
    val params.GENOME_PREFIX

    output:
    file 'out.gg' into gg_file_ref
    file 'shortname' into shortname_ref
    file 'mapped.fasta' into mapped_fasta_ref

    """
    truncate_header.lua < pepfile.fas > pepfile.trunc
    ln -s pepfile.trunc mapped.fasta
    make_gg_line.lua ${params.GENOME_PREFIX} mapped.fasta  > out.gg
    echo "${params.GENOME_PREFIX}" > shortname
    """
}

full_gg = gg_file.mix(gg_file_ref).collectFile()
full_shortnames = shortname.mix(shortname_ref).collectFile()
full_mapped_fasta = mapped_fasta.mix(mapped_fasta_ref).collectFile()

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
    makeblastdb -dbtype prot -in mapped.fasta
    """
}

proteins_orthomcl_blast_chunk = full_mapped_fasta_for_query.splitFasta( by: 50)
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
   # blastp -word_size 6 -evalue 1e-5 -db mapped.fasta -outfmt 6 \
   #  -query mapped_chunk.fasta > blastout
    blastall -p blastp -W 4 -e 0.00001 -F T -d mapped.fasta -m 8 \
      -i mapped_chunk.fasta > blastout
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

orthomcl_cluster_out_annot = Channel.create()
result_ortho = Channel.create()
orthomcl_cluster_out.into(orthomcl_cluster_out_annot, result_ortho)

process annotate_orthologs {
    cache 'deep'

    input:
    file 'orthomcl_out' from orthomcl_cluster_out_annot
   //  file 'mapfile' from full_mapfile
    file 'input.gff3' from genemodels_for_omcl_annot
    file 'gff_ref.gff3' from omcl_gfffile
    file 'gaf_ref.gaf' from omcl_gaffile

    output:
    file 'with_func.gff3' into gff3_with_ortho_transferred
    file 'orthomcl.gaf' into gaf_with_ortho_transferred

    """
    # annotate GFF with ortholog clusters and members
    map_clusters_gff.lua input.gff3 orthomcl_out > with_clusters.gff3

    # transfer functional annotation from orthologs
    transfer_annotations_from_gff.lua with_clusters.gff3 \
        gff_ref.gff3 > with_func.gff3

    # transfer GOs from orthologs GAF
    transfer_annotations_from_gaf.lua with_func.gff3 \
        gaf_ref.gaf ${params.DB_ID} \
        ${params.TAXON_ID} > orthomcl.gaf
    """
}

// PFAM DOMAIN ANNOTATION
// ======================

proteins_pfam_chunk = proteins_pfam.splitFasta( by: 30)
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
    pfam_to_gff3.lua ${PFAM2GO} < pfam.out \
      | gt gff3 -sort -tidy -retainids > pfam.gff3
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

// MAKE DISTRIBUTION
// =================

process make_distribution_gff {
    input:
    file 'pseudo.in.gff3' from gff3_with_pfam
    file 'pseudo.pseudochr.agp' from pseudochr_agp_make_dist
    file 'pseudo.scafs.agp' from scaffolds_agp_make_dist
    file 'pseudochr.fasta' from pseudochr_seq_make_dist_1
    file 'scafs.fasta' from scaffolds_seq_make_dist_1

    output:
    set file('pseudo.out.gff3'),
        file('scaffold.out.gff3') into result_gff3
    set file('pseudo.pseudochr.agp'),
        file('pseudo.scafs.agp') into result_agp

    """
    split_gff_by_agp.lua pseudo.in.gff3 pseudo.pseudochr.agp \
      pseudochr.fasta scafs.fasta > 1
    gt gff3 -sort -retainids -tidy 1 > scaffold.out.gff3
    cp pseudo.in.gff3 pseudo.out.gff3
    """
}

process make_distribution_gaf {
    input:
    file 'in.gaf' from gaf_with_pfam

    output:
    file 'out.gaf' into result_gaf

    """
    echo '!gaf-version: 1.0' > out.gaf
    sort -k2,2 in.gaf >> out.gaf
    """
}

process make_distribution_seqs {
    input:
    file 'pseudochr.in' from pseudochr_seq_make_dist_2
    file 'scafs.in' from scaffolds_seq_make_dist_2

    output:
    set file('pseudochr.fasta.gz'),
        file('scafs.fasta.gz') into result_seq

    """
    cp pseudochr.in ./pseudochr.fasta
    gzip pseudochr.fasta
    cp scafs.in ./scafs.fasta
    gzip scafs.fasta
    """
}

stats_inseq = Channel.create()
circos_inseq = Channel.create()
report_inseq = Channel.create()
embl_inseq = Channel.create()
out_seq = Channel.create()
result_seq.into(stats_inseq, circos_inseq, report_inseq, out_seq, embl_inseq)

stats_gff3 = Channel.create()
circos_gff3 = Channel.create()
report_gff3 = Channel.create()
embl_gff3 = Channel.create()
out_gff3 = Channel.create()
refcomp_gff3_in = Channel.create()
genelist_gff3_in = Channel.create()
prot_fasta_annot_gff3 = Channel.create()
result_gff3.into(stats_gff3, circos_gff3, report_gff3, out_gff3, embl_gff3,
                 refcomp_gff3_in, genelist_gff3_in, prot_fasta_annot_gff3)

// GENOME STATS GENERATION
// =======================

process make_genome_stats {
    input:
    set file('pseudo.fasta.gz'), file('scaf.fasta.gz') from stats_inseq
    set file('pseudo.gff3'), file('scaf.gff3') from stats_gff3

    output:
    file 'stats.txt' into stats_output

    """
    gunzip -f pseudo.fasta.gz
    genome_stats.lua pseudo.gff3 pseudo.fasta > stats.txt && rm -f pseudo.fasta
    """
}

// CIRCOS PLOTS
// ============

if (params.do_contiguation && params.do_circos) {
    process blast_for_circos {
        input:
        set file('pseudo.fasta.gz'), file('scaf.fasta.gz') from circos_inseq
        file 'refseq.fasta' from ref_seq

        output:
        file 'blastout.txt' into circos_blastout

        """
        gunzip -f pseudo.fasta.gz
        makeblastdb -dbtype nucl -in refseq.fasta
        blastn -db refseq.fasta -query pseudo.fasta -evalue 1e-6 -outfmt 6 \
          -out blastout.txt
        """
    }

    process make_circos_inputs {
        input:
        set file('pseudo.gff3'), file('scaf.gff3') from circos_gff3
        file 'refannot.gff3' from ref_annot
        file 'blast.in' from circos_blastout
        val params.CHR_PATTERN
        val params.ABACAS_BIN_CHR

        output:
        file 'links.txt' into circos_input_links
        file 'karyotype.txt' into circos_input_karyotype
        file 'chromosomes.txt' into circos_input_chromosomes
        file 'genes.txt' into circos_input_genes
        file 'gaps.txt' into circos_input_gaps
        file 'bin.txt' into bin_target_mapping

        """
        prepare_circos_inputs.lua refannot.gff3 pseudo.gff3 blast.in . \
           "${params.CHR_PATTERN}" "${params.ABACAS_BIN_CHR}"
        """
    }

    circos_chromosomes = ref_target_mapping.splitCsv(sep: "\t")
    circos_conffile = file(params.CIRCOS_CONFIG_FILE)
    process circos_run_chrs {
        tag { chromosome[0] }

        input:
        file 'links.txt' from circos_input_links.first()
        file 'karyotype.txt' from circos_input_karyotype.first()
        file 'genes.txt' from circos_input_genes.first()
        file 'gaps.txt' from circos_input_gaps.first()
        val circos_conffile
        val chromosome from circos_chromosomes

        output:
        set file('image.png'), val(chromosome) into circos_output

        """
        circos  -conf ${circos_conffile} -param image/file=image.png  \
                -param chromosomes='${chromosome[1]};${chromosome[2]}' \
                -param chromosomes_reverse=${chromosome[1]}
        """
    }

    circos_binmap = bin_target_mapping.splitCsv(sep: "\t")
    process circos_run_bin {
        tag { "bin" }

        input:
        file 'links.txt' from circos_input_links.first()
        file 'karyotype.txt' from circos_input_karyotype.first()
        file 'genes.txt' from circos_input_genes.first()
        file 'gaps.txt' from circos_input_gaps.first()
        val circos_conffile
        val binmap from circos_binmap

        output:
        set file('bin.png'), val('bin') into circos_output

        script:
        if (binmap[1] && binmap[2])
            """
            circos  -conf ${circos_conffile} -param image/file=bin.png  \
                    -param chromosomes='${binmap[1]};${binmap[2]}' \
                    -param chromosomes_reverse=${binmap[1]}
            """
    }

    circos_output.subscribe {
        println it[0]
        if (params.dist_dir) {
          it[0].copyTo(params.dist_dir + "/chr" + it[1][0] + ".png")
        }
    }
}

// EMBL OUTPUT
// ===========

if (params.make_embl) {
    process merge_gff3_for_gff3toembl {
        input:
        set file('pseudo.fasta.gz'), file('scaf.fasta.gz') from embl_inseq
        set file('pseudo.gff3'), file('scaf.gff3') from embl_gff3

        output:
        file 'full.gff3' into embl_full_gff

        """
        gt inlineseq_add -seqfile pseudo.fasta.gz -matchdescstart -force \
          -o full.gff3 pseudo.gff3
        """
    }

    process make_embl {
        input:
        file 'embl_in.gff3' from embl_full_gff

        output:
        file 'outfile.embl' into embl_out

        """
        gff3_to_embl --authors '${params.EMBL_AUTHORS}' \
          --title '${params.EMBL_TITLE}' \
          --publication '${params.EMBL_PUBLICATION}' \
          --genome_type '${params.EMBL_GENOME_TYPE}' \
          --classification '${params.EMBL_CLASSIFICATION}' \
          --output_filename outfile.embl \
          '${params.EMBL_ORGANISM}' '${params.TAXON_ID}' \
            '${params.EMBL_PROJ_ACCESSION}' '${params.EMBL_DESCRIPTION}' \
            embl_in.gff3
        """
    }

    embl_out.subscribe {
        println it
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }
}

// REFERENCE COMPARISON
// ====================

if (params.use_reference) {
    process reference_compare {
        input:
        file 'in.protein.fasta' from refcomp_protein_in
        set file('pseudo.in.annotation.gff3'), file('scafs.in.annotation.gff3') from refcomp_gff3_in
        val params.GENOME_PREFIX
        val params.ref_dir

        output:
        file 'tree_selection.fasta' into tree_fasta
        file 'tree_selection.genes' into tree_genes
        file 'in.protein.fasta' into refcomp_protein_out

        """
        stream_new_against_core.lua pseudo.in.annotation.gff3 in.protein.fasta \
          ${params.ref_dir} '' ${params.GENOME_PREFIX}
        """
    }

    // TREE CALCULATION
    // ================

    process make_tree {
        input:
        file 'tree_selection.fasta' from tree_fasta

        output:
        file "tree.out" into tree_out
        file "tree.aln" into tree_aln

        """
        mafft --auto --anysymbol --parttree --quiet tree_selection.fasta > tree.aln
        FastTree tree.aln > tree.out
        """
    }

    tree_out.subscribe {
        println it
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }

    tree_genes.subscribe {
        println it
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }

    tree_aln.subscribe {
        println it
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }
}

// REPORT CREATION
// ===============

specfile = file(params.SPECFILE)
process make_report {
    input:
    set file('pseudo.fasta.gz'), file('scaf.fasta.gz') from report_inseq
    set file('pseudo.gff3'), file('scaf.gff3') from report_gff3
    val specfile

    output:
    set file('pseudo.report.html'), file('scaf.report.html') into report_output

    """
    gt speck -specfile ${specfile} -matchdescstart -seqfile pseudo.fasta.gz \
      -provideindex -typecheck so -output html pseudo.gff3 > pseudo.report.html
    gt speck -specfile ${specfile} -matchdescstart -seqfile scaf.fasta.gz \
      -provideindex -typecheck so -output html scaf.gff3 > scaf.report.html
    """
}

// SIMPLE GENELIST
// ===============

process make_genelist {
    input:
    set file('pseudo.gff3'), file('scaf.gff3') from genelist_gff3_in

    output:
    file 'genelist.csv' into genelist_csv_out

    """
    genes_gff3_to_csv.lua pseudo.gff3 > genelist.csv
    """
}

genelist_csv_out.subscribe {
    println it
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

// PROTEIN SEQUENCE IMPROVEMENT
// ============================

process add_products_to_protein_fasta {
    input:
    file 'in.fasta' from proteins_output
    set file('pseudo.gff3'), file('scaf.gff3') from prot_fasta_annot_gff3

    output:
    file 'proteins.fasta' into result_protein

    """
    add_products_to_fasta.lua pseudo.gff3 in.fasta > proteins.fasta
    """
}

// OUTPUT
// ======

out_gff3.subscribe {
    println it
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}

out_seq.subscribe {
    println it
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}

result_agp.subscribe {
    println it
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}

result_gaf.collectFile().subscribe {
    println it
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

result_ortho.collectFile().subscribe {
    println it
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

result_protein.subscribe {
    println it
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

stats_output.subscribe {
    println it
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

report_output.subscribe {
    println it
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}