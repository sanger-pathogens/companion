#!/usr/bin/env nextflow

VERSION = "1.0.2"

/*
    Author: Sascha Steinbiss <ss34@sanger.ac.uk>
    Copyright (c) 2014-2016 Genome Research Ltd

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

genome_file = file(params.inseq)
ref_annot = file(params.ref_dir + "/" + params.ref_species + "/annotation.gff3")
ref_seq = file(params.ref_dir + "/" + params.ref_species + "/genome.fasta")
ref_dir = file(params.ref_dir)
go_obo = file(params.GO_OBO)
ncrna_models = file(params.NCRNA_MODELS)
extrinsic_cfg = file(params.AUGUSTUS_EXTRINSIC_CFG)
omcl_gfffile = file(params.ref_dir + "/" + params.ref_species + "/annotation.gff3")
omcl_gaffile = file(params.ref_dir + "/" + params.ref_species + "/go.gaf")
omcl_pepfile = file(params.ref_dir + "/" + params.ref_species + "/proteins.fasta")
augustus_modeldir = file(params.ref_dir + "/" + params.ref_species)

log.info ""
log.info "C O M P A N I O N  ~  version " + VERSION
log.info "query               : ${params.inseq}"
log.info "reference           : ${params.ref_species}"
log.info "reference directory : ${params.ref_dir}"
if (params.dist_dir) {
    distDir = new File(params.dist_dir)
    if(distDir.exists() && distDir.isFile()) {
        exit 1, "cannot create output path: ${params.dist_dir} -- already exists as file"
    }
    if(!distDir.exists() && !distDir.mkdirs()) {
        exit 1, "cannot create output path: ${params.dist_dir} -- check file system permissions"
    }
    if(!distDir.isDirectory()) {
        exit 1, "cannot prepare output path: ${params.dist_dir} -- aborting"
    }
    log.info "output directory    : ${params.dist_dir}"
}
log.info ""

// INPUT SANITIZATION
// ==================

if (params.truncate_input_headers) {
    process truncate_input_headers {

        input:
        file genome_file

        output:
        file 'truncated.fasta' into truncated_genome_file

        """
        truncate_header.lua < ${genome_file} > truncated.fasta
        """
    }
} else {
    genome_file.into { truncated_genome_file }
}

process sanitize_input {
    input:
    file 'truncated.fasta' from truncated_genome_file

    output:
    file 'sanitized.fasta' into sanitized_genome_file

    """
    sed 's/['\\''+&= ]/_/g' truncated.fasta | trim_wildcards.lua > sanitized.fasta
    """
}

// PSEUDOCHROMOSOME CONTIGUATION
// =============================

if (params.do_contiguation) {
    ref_chr = file(params.ref_dir + "/" + params.ref_species + "/chromosomes.fasta")
    process contiguate_pseudochromosomes {
        afterScript 'rm -rf Ref.* Res.*'

        input:
        file sanitized_genome_file
        file ref_chr
        file ref_dir

        output:
        file 'pseudo.pseudochr.fasta' into pseudochr_seq
        file 'pseudo.pseudochr.agp' into pseudochr_agp
        file 'pseudo.scafs.fasta' into scaffolds_seq
        file 'pseudo.scafs.agp' into scaffolds_agp
        file 'pseudo.contigs.fasta' into contigs_seq
        file 'ref_target_mapping.txt' into ref_target_mapping

        """
        abacas2.nonparallel.sh \
          "${ref_chr}" "${sanitized_genome_file}" "${params.ABACAS_MATCH_SIZE}" \
          "${params.ABACAS_MATCH_SIM}" 0 3000
        abacas_combine.lua . pseudo "${ref_dir}" "${params.ref_species}" \
          "${params.GENOME_PREFIX}" "${params.ABACAS_BIN_CHR}" \
          "${params.GENOME_PREFIX}"
        """
    }
} else {
    process prepare_noncontiguated_input {
        input:
        file sanitized_genome_file

        output:
        file 'pseudo.pseudochr.fasta' into pseudochr_seq
        file 'pseudo.pseudochr.agp' into pseudochr_agp
        file 'pseudo.scafs.fasta' into scaffolds_seq
        file 'pseudo.scafs.agp' into scaffolds_agp
        file 'pseudo.contigs.fasta' into contigs_seq

        """
        no_abacas_prepare.lua ${sanitized_genome_file} pseudo
        """
    }
}

// fork important output streams
pseudochr_seq.into{ pseudochr_seq_tRNA
					pseudochr_seq_ncRNA
					pseudochr_seq_ratt
					pseudochr_seq_augustus
					pseudochr_seq_augustus_ctg
                    pseudochr_seq_snap
                    pseudochr_seq_make_gaps
                    pseudochr_seq_make_dist_1
                    pseudochr_seq_make_dist_2
                    pseudochr_seq_splitsplice
                    pseudochr_seq_pseudogene
                    pseudochr_seq_integrate
                    pseudochr_seq_tmhmm
                    pseudochr_seq_orthomcl
                    pseudochr_seq_exonerate }

scaffolds_seq.into{ scaffolds_seq_augustus
				    scaffolds_seq_make_gaps
                    scaffolds_seq_make_dist_1
                    scaffolds_seq_make_dist_2 }


scaffolds_agp.into{ scaffolds_agp_augustus
                    scaffolds_agp_rnaseq
					scaffolds_agp_make_gaps
                    scaffolds_agp_make_dist }

pseudochr_agp.into{ pseudochr_agp_augustus
                    pseudochr_agp_rnaseq
				    pseudochr_agp_make_gaps
                    pseudochr_agp_make_dist }

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

process press_ncRNA_cms {
    input:
    val ncrna_models

    output:
    file 'models.cm*' into ncrna_cmindex

    """
    cp ${ncrna_models} ./models.cm
    cmpress -F models.cm
    """
}

pseudochr_seq_ncRNA.splitFasta( by: 5, file: true).set { ncrna_genome_chunk }

process predict_ncRNA {
    input:
    file 'chunk' from ncrna_genome_chunk
    file ncrna_models from ncrna_cmindex.first()

    output:
    file 'cm_out' into cmtblouts

    """
    cmsearch --cpu 1 --tblout cm_out --cut_ga models.cm chunk
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
    process make_exonerate_index {
        input:
        file 'genome.fasta' from pseudochr_seq_exonerate.first()

        output:
        file 'index.esi' into exn_index_esi
        file 'index.esd' into exn_index_esd
        file 'genome.fasta' into exn_index_fasta

        """
        fasta2esd --softmask no genome.fasta index.esd
        esd2esi index.esd index.esi --translate yes
        """

    }

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

    exn_prot_chunk = ref_pep.splitFasta( by: 200, file: true)
    process run_exonerate {
        cache 'deep'
        // this process can fail for rogue exonerate processes
        errorStrategy 'ignore'
        time '3h'

        input:
        file 'index.esi' from exn_index_esi.first()
        file 'index.esd' from exn_index_esd.first()
        file 'genome.fasta' from exn_index_fasta.first()
        file 'prot.fasta' from exn_prot_chunk

        output:
        file 'exn_out' into exn_results

        """
        get_unused_port.sh > port
        reaper.sh exonerate-server --port `cat port` --input index.esi &
        sleep 10
        exonerate -E false --model p2g --showvulgar no --showalignment no \
          --showquerygff no --showtargetgff yes --percent 80 --geneseed 250 \
          --ryo \"AveragePercentIdentity: %pi\n\" prot.fasta \
          localhost:`cat port` > exn_out
        kill `cat exonerate-server.pid`
        rm -f exonerate-server.pid
        """
    }

    process exonerate_make_hints {
        cache 'deep'

        input:
        file 'exnout' from exn_results.collectFile()

        output:
        file 'augustus.hints' into exn_hints

        script:
        """
        exonerate2hints.pl \
          --source=P --maxintronlen=${params.AUGUSTUS_HINTS_MAXINTRONLEN} \
          --in=${exnout} \
          --out=augustus.hints
        """
    }
} else {
    process exonerate_empty_hints {
        output:
        file 'augustus.hints' into exn_hints

        """
        touch augustus.hints
        """
    }
}

// RATT
// ====

if (params.run_ratt) {
    process ratt_make_ref_embl {
        input:
        file ref_annot
        file ref_seq
        val go_obo

        output:
        file '*.embl' into ref_embl

        """
        gff3_to_embl.lua -o ${ref_annot} ${go_obo} Foo ${ref_seq}
        """
    }

    process run_ratt {
        afterScript 'rm -rf Reference* Sequences Query query.*'

        input:
        file 'in*.embl' from ref_embl
        file 'pseudo.pseudochr.fasta' from pseudochr_seq_ratt

        output:
        file 'Out*.final.embl' into ratt_result
        file 'Out*.Report.txt' into ratt_reports

        """
        touch Out.0.Report.txt
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
        echo '##gff-version 3' > ratt.gff3
        ratt_embl_to_gff3.lua in*.embl | \
          gt gff3 -sort -retainids -tidy > \
          ratt.tmp.gff3
        if [ -s ratt.tmp.gff3 ]; then
          ratt_remove_problematic.lua ratt.tmp.gff3 in*report | \
          gt gff3 -sort -retainids -tidy > ratt.gff3;
        fi
        """
    }
} else {
    process ratt_empty_models {
    output:
    file 'result.gff3' into ratt_gff3

    """
    echo '##gff-version 3' > result.gff3
    """
  }
}

// TRANSCRIPT EVIDENCE PREPARATION
// ===============================

if (params.TRANSCRIPT_FILE) {
  transcript_evidence = file(params.TRANSCRIPT_FILE)

  process transform_input_gtf {
    input:
    file 'transcripts.gtf' from transcript_evidence
    file 'pseudochr.agp' from pseudochr_agp_rnaseq

    output:
    file 'transcripts_transformed.gtf' into transcript_evidence_transformed

    """
    LC_ALL='C' sort -k 1,1 -k 4,4n transcripts.gtf | uniq > 1
    transform_gtf_with_agp.lua 1 pseudochr.agp > transcripts_transformed.gtf
    """
  }

  process prepare_transcript_hints {
    input:
    file 'transcripts.gtf' from transcript_evidence_transformed

    output:
    file 'transcripts.hints' into trans_hints

    """
    LC_ALL='C' sort -k 1,1 -k 4,4n transcripts.gtf | uniq | \
      cufflinks_to_hints.lua > transcripts.hints
    """
  }
} else {
  process transcript_empty_hints {
    output:
    file 'transcripts.hints' into trans_hints

    """
    touch transcripts.hints
    """
  }
}

// HINTS PREPARATION
// =================

process merge_hints {
    cache 'deep'

    input:
    file 'hints.exon.txt' from exn_hints
    file 'hints.trans.txt' from trans_hints

    output:
    set stdout, file('hints.txt') into all_hints

    """
    touch hints.txt
    cat hints.exon.txt hints.trans.txt > hints.concatenated.txt
    if [ -s hints.concatenated.txt ] ; then
      mv hints.concatenated.txt hints.txt;
      echo -n '--alternatives-from-evidence=false --hintsfile=augustus.hints';
    fi
    """
}

// AUGUSTUS
// ========

process run_augustus_pseudo {
    cache 'deep'

    input:
    set val(hintsline), file('augustus.hints') from all_hints
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_augustus
    val extrinsic_cfg
    file augustus_modeldir

    output:
    file 'augustus.gff3' into augustus_pseudo_gff3

    """
    echo "##gff-version 3\n" > augustus.full.tmp.2;
    AUGUSTUS_CONFIG_PATH=${augustus_modeldir} \
        augustus \
            --species=augustus_species \
            --stopCodonExcludedFromCDS=false \
            --protein=off --codingseq=off --strand=both \
            --genemodel=${params.AUGUSTUS_GENEMODEL} --gff3=on \
            ${hintsline} \
            --noInFrameStop=true \
            --extrinsicCfgFile=${extrinsic_cfg} \
            pseudo.pseudochr.fasta > augustus.full.tmp
    augustus_to_gff3.lua < augustus.full.tmp \
        | gt gff3 -sort -tidy -retainids > 1
    if [ -s 1 ]; then
        gt select -mingenescore ${params.AUGUSTUS_SCORE_THRESHOLD} 1 \
        > augustus.full.tmp.2;
    fi
    augustus_mark_partial.lua augustus.full.tmp.2 > augustus.gff3
    """
}

process run_augustus_contigs {
    input:
    file 'pseudo.contigs.fasta' from contigs_seq
    file 'pseudo.scaffolds.agp' from scaffolds_agp_augustus
    file 'pseudo.scaffolds.fasta' from scaffolds_seq_augustus
    file 'pseudo.pseudochr.agp' from pseudochr_agp_augustus
    file 'pseudo.pseudochr.fasta' from pseudochr_seq_augustus_ctg
    file augustus_modeldir

    output:
    file 'augustus.scaf.pseudo.mapped.gff3' into augustus_ctg_gff3

    """
    echo "##gff-version 3\n" > augustus.ctg.tmp.2;
    AUGUSTUS_CONFIG_PATH=${augustus_modeldir} \
        augustus --species=augustus_species \
            --stopCodonExcludedFromCDS=false \
            --protein=off --codingseq=off --strand=both --genemodel=partial \
            --gff3=on \
            --noInFrameStop=true \
            pseudo.contigs.fasta > augustus.ctg.tmp
    augustus_to_gff3.lua < augustus.ctg.tmp \
        | gt gff3 -sort -tidy -retainids > 1
    if [ -s 1 ]; then
        gt select -mingenescore ${params.AUGUSTUS_SCORE_THRESHOLD} 1 \
        > augustus.ctg.tmp.2;
    fi
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

    if [ ! -s augustus.scaf.pseudo.mapped.gff3 ]; then
      echo '##gff-version 3' > augustus.scaf.pseudo.mapped.gff3
    fi
    """
}

if (params.run_snap) {
    snap_model = file(params.ref_dir + "/" + params.ref_species + "/snap.hmm")
    process run_snap {
        input:
        file 'pseudo.pseudochr.fasta' from pseudochr_seq_snap
        file 'snap.hmm' from snap_model

        output:
        file 'snap.gff3' into snap_gff3

        """
        echo '##gff-version 3' > snap.gff3
        snap -gff -quiet snap.hmm pseudo.pseudochr.fasta > snap.tmp
        snap_gff_to_gff3.lua snap.tmp > snap.tmp.2
        if [ -s 1 ]; then
            gt gff3 -sort -tidy -retainids snap.tmp.2 > snap.gff3;
        fi
        """
    }
} else {
    process make_empty_snap {
        output:
        file 'snap.gff3' into snap_gff3

        """
        echo '##gff-version 3' > snap.gff3
        """
    }
}

process merge_genemodels {
    cache 'deep'

    input:
    file 'augustus.full.gff3' from augustus_pseudo_gff3
    file 'augustus.ctg.gff3' from augustus_ctg_gff3
    file 'snap.full.gff3' from snap_gff3
    file 'ratt.full.gff3' from ratt_gff3

    output:
    file 'merged.gff3' into merged_gff3

    """
    unset GT_RETAINIDS && \
    gt gff3 -fixregionboundaries -retainids no -sort -tidy \
        augustus.full.gff3 augustus.ctg.gff3 snap.full.gff3 ratt.full.gff3 \
        > merged.pre.gff3 && \
    export GT_RETAINIDS=yes
    if [ ! -s merged.pre.gff3 ]; then
        echo '##gff-version 3' > merged.pre.gff3
    fi

    # avoid huge gene clusters
    gt select -maxgenelength ${params.MAX_GENE_LENGTH} merged.pre.gff3 > merged.gff3
    """
}

process integrate_genemodels {
    cache 'deep'
    afterScript 'rm -rf sequence.fasta.*'

    input:
    file 'merged.gff3' from merged_gff3
    file 'sequence.fasta' from pseudochr_seq_integrate

    output:
    file 'integrated.gff3' into integrated_gff3

    script:
    if (params.WEIGHT_FILE.length() > 0)
        """
        integrate_gene_calls.lua -w ${params.WEIGHT_FILE} -s sequence.fasta < merged.gff3 | \
            gt gff3 -sort -tidy -retainids > integrated.gff3
        """
    else
        """
        integrate_gene_calls.lua -s sequence.fasta < merged.gff3 | \
            gt gff3 -sort -tidy -retainids > integrated.gff3
        """
}

if (params.fix_polycistrons) {
    process fix_polycistrons {
        input:
        file 'integrated.gff3' from integrated_gff3

        output:
        file 'integrated.fixed.sorted.gff3' into integrated_gff3_processed

        """
        if [ ! -s integrated.gff3 ]; then
          echo '##gff-version 3' > integrated.gff3
        fi

        fix_polycistrons.lua integrated.gff3 > integrated.fixed.gff3

        # make sure final output is sorted
        gt gff3 -sort -tidy -retainids \
          integrated.fixed.gff3 > integrated.fixed.sorted.gff3
        """
    }
} else {
    integrated_gff3.into { integrated_gff3_processed }
}

process remove_exons {
    input:
    file 'integrated.gff3' from integrated_gff3_processed

    output:
    file 'integrated_clean.gff3' into integrated_gff3_clean

    """
    if [ ! -s integrated.gff3 ]; then
      echo '##gff-version 3' > integrated.gff3
    fi
    remove_exons.lua integrated.gff3 > integrated_clean.gff3
    """
}

if (params.do_pseudo) {
    process pseudogene_indexing {
        input:
        file 'ref.peps.fasta' from omcl_pepfile

        output:
        file 'prot_index*' into pseudochr_last_index

        """
        tantan -p -r0.02 ref.peps.fasta | lastdb -p -c prot_index
        """
    }

    pseudochr_seq_pseudogene.into{ pseudochr_seq_pseudogene_align; pseudochr_seq_pseudogene_calling }
    pseudochr_seq_pseudogene_align.splitFasta( by: 3, file: true).set { pseudogene_align_chunk }

    process pseudogene_last {
        input:
        file 'chunk.fasta' from pseudogene_align_chunk
        file prot_index from pseudochr_last_index.first()

        output:
        file 'last.out' into pseudochr_last_out

        """
        lastal -R01 -pBL80 -F15 -e400 -m10 -f0 prot_index chunk.fasta > last.out
        """
    }

    process pseudogene_calling {
        cache 'deep'

        input:
        file 'pseudochr.fasta' from pseudochr_seq_pseudogene_calling
        file 'genes.gff3' from integrated_gff3_clean
        file 'last.out' from pseudochr_last_out.collectFile()

        output:
        file 'genes_and_pseudo.gff3' into gff3_with_pseudogenes

        """
        # reconstruct frameshifted candidates from output
        pseudo_merge_last.lua last.out pseudochr.fasta > last_gff.gff3

        # merge with gene models
        gt gff3 -sort -tidy -retainids last_gff.gff3 genes.gff3 > last_and_genes.gff3
        if [ ! -s last_and_genes.gff3 ]; then
          echo '##gff-version 3' > last_and_genes.gff3
        fi
        pseudo_merge_with_genes.lua last_and_genes.gff3 pseudochr.fasta > out_tmp.gff3
        gt gff3 -sort -retainids -tidy out_tmp.gff3 > genes_and_pseudo.gff3
        """
    }
} else {
    gff3_with_pseudogenes = integrated_gff3_clean
}

// MERGE ALL GENES TO FINAL SET AND CLEANUP
// ========================================

process merge_structural {
    cache 'deep'

    input:
    file 'ncrna.gff3' from ncrnafile
    file 'trna.gff3' from trnas
    file 'integrated.gff3' from gff3_with_pseudogenes

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
      gt gff3 -fixregionboundaries -sort -tidy -retainids > contigs.gff3

    transform_gff_with_agp.lua contigs.gff3 \
      pseudo.pseudochr.agp pseudo.scaffolds.fasta \
      pseudo.pseudochr.fasta "between scaffolds" | \
      gt gff3 -fixregionboundaries -sort -tidy -retainids > contigs2.gff3

    if [ ! -s merged_in.gff3 ]; then
      echo '##gff-version 3' > merged_in.gff3
    fi

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
    splice_genes_at_gaps_new.lua tmp3 | \
      gt gff3 -sort -retainids -tidy | \
      gt select -rule_files ${FILTER_SHORT_PARTIALS_RULE} > tmp4

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

    script:
    if (params.alphanumeric_ids)
        """
        create_polypeptides.lua -r input.gff3 ${params.GENOME_PREFIX} \
            "${params.CHR_PATTERN}" | gt gff3 -sort -retainids -tidy > output.gff3
        """
    else
        """
        create_polypeptides.lua input.gff3 ${params.GENOME_PREFIX} \
            "${params.CHR_PATTERN}" | gt gff3 -sort -retainids -tidy > output.gff3
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

genemodels_with_polypeptides_gff3.into{ genemodels_for_omcl_proteins; genemodels_for_omcl_annot }

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

proteins_target.into{ proteins_orthomcl
					  proteins_pfam
					  refcomp_protein_in
                      proteins_output }

process make_ref_input_for_orthomcl {
    input:
    file omcl_pepfile

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

gg_file.mix(gg_file_ref).collectFile().set { full_gg }
shortname.mix(shortname_ref).collectFile().set { full_shortnames }
mapped_fasta.mix(mapped_fasta_ref).collectFile().set { full_mapped_fasta }

full_mapped_fasta.into{ full_mapped_fasta_for_index; full_mapped_fasta_for_query }

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

proteins_orthomcl_blast_chunk = full_mapped_fasta_for_query.splitFasta( by: 50, file: true)
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
    blastall -p blastp -W 4 -F 'm S' -v 100000 -b 100000 -d mapped.fasta -m 8 \
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

orthomcl_cluster_out.into{ orthomcl_cluster_out_annot; result_ortho }

process annotate_orthologs {
    cache 'deep'

    input:
    file 'orthomcl_out' from orthomcl_cluster_out_annot
    file 'input.gff3' from genemodels_for_omcl_annot
    file 'gff_ref.gff3' from omcl_gfffile
    file 'gaf_ref.gaf' from omcl_gaffile

    output:
    file 'with_func.gff3' into gff3_with_ortho_transferred
    file 'orthomcl.gaf' into gaf_with_ortho_transferred

    """
    # ensure sorting of input file
    gt gff3 -sort -retainids -tidy input.gff3 > input.gff3.sorted

    # annotate GFF with ortholog clusters and members
    map_clusters_gff.lua input.gff3.sorted orthomcl_out > with_clusters.gff3

    # transfer functional annotation from orthologs
    transfer_annotations_from_gff.lua with_clusters.gff3 \
        gff_ref.gff3 > with_func.gff3

    # transfer GOs from orthologs GAF
    transfer_annotations_from_gaf.lua with_func.gff3 \
        gaf_ref.gaf ${params.DB_ID} \
        ${params.TAXON_ID} > orthomcl.gaf

    # clean up
    rm -f input.gff3 input.gff3.sorted with_clusters.gff3
    """
}

// PFAM DOMAIN ANNOTATION
// ======================

proteins_pfam_chunk = proteins_pfam.splitFasta( by: 30, file: true)
process run_pfam {
    cache 'deep'

    input:
    file 'proteins.fas' from proteins_pfam_chunk

    output:
    file 'pfamout' into pfam_output

    """
    hmmscan --domtblout pfamout --cut_ga --noali ${PFAM} proteins.fas
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
    if [ ! -s pfam.gff3 ]; then
        echo '##gff-version 3' > pfam.gff3
    fi
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
    cp pseudochr.in ./pseudochr.fasta && gzip -9 pseudochr.fasta
    cp scafs.in ./scafs.fasta && gzip -9 scafs.fasta
    """
}

result_seq.into{ stats_inseq
				 circos_inseq
				 report_inseq
				 out_seq
				 embl_inseq }


result_gff3.into{ stats_gff3
				  circos_gff3
				  report_gff3
				  out_gff3
				  embl_gff3
                  refcomp_gff3_in
                  genelist_gff3_in
                  prot_fasta_annot_gff3 }

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

// REFERENCE COMPARISON
// ====================

if (params.use_reference) {
    process reference_compare {
        input:
        file 'in.protein.fasta' from refcomp_protein_in
        set file('pseudo.in.annotation.gff3'), file('scafs.in.annotation.gff3') from refcomp_gff3_in
        file ref_dir

        output:
        file 'tree_selection.fasta' into tree_fasta
        file 'tree_selection.genes' into tree_genes
        file 'in.protein.fasta' into refcomp_protein_out
        file 'core_comp_circos.txt' into core_comp_circos
        file 'core_comparison.txt' into core_comp

        """
        stream_new_against_core.lua pseudo.in.annotation.gff3 in.protein.fasta \
          ${ref_dir} ${params.GENOME_PREFIX} ${params.ref_species}
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
        touch tree.out
        touch tree.aln
        if [ -s tree_selection.fasta ] ; then
          mafft --auto --anysymbol --retree 1 --parttree --quiet tree_selection.fasta > tree.aln;
          FastTree tree.aln > tree.out;
        fi
        """
    }

    tree_out.subscribe {
        if (params.print_paths) {
          println it
        }
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }

    tree_genes.subscribe {
        if (params.print_paths) {
          println it
        }
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }

    tree_aln.subscribe {
        if (params.print_paths) {
          println it
        }
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }

    core_comp.subscribe {
        if (params.print_paths) {
          println it
        }
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }
} else {
    process make_empty_circos_clusters {
        output:
        file 'core_comp_circos.txt' into core_comp_circos

        """
        touch core_comp_circos.txt
        """
    }
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
        file ref_dir

        output:
        file 'links.txt' into circos_input_links
        file 'karyotype.txt' into circos_input_karyotype
        file 'chromosomes.txt' into circos_input_chromosomes
        file 'genes.txt' into circos_input_genes
        file 'gaps.txt' into circos_input_gaps
        file 'bin.txt' into bin_target_mapping

        """
        prepare_circos_inputs.lua refannot.gff3 pseudo.gff3 blast.in . \
           "${params.CHR_PATTERN}" "${params.ABACAS_BIN_CHR}" "${ref_dir}" \
           "${params.ref_species}"
        """
    }

    circos_input_links.into{ circos_input_links_chr; circos_input_links_bin }
    core_comp_circos.into{ core_comp_circos_chr; core_comp_circos_bin }

    circos_input_karyotype.into{ circos_input_karyotype_chr; circos_input_karyotype_bin }

    circos_input_chromosomes.into{ circos_input_chromosomes_chr; circos_input_chromosomes_bin }

    circos_input_genes.into{ circos_input_genes_chr; circos_input_genes_bin }

    circos_input_gaps.into{ circos_input_gaps_chr; circos_input_gaps_bin }

    ref_target_mapping.splitCsv(sep: "\t").set { circos_chromosomes }
    circos_conffile = file(params.CIRCOS_CONFIG_FILE)

    process circos_run_chrs {
        tag { chromosome[0] }

        input:
        file 'links.txt' from circos_input_links_chr.first()
        file 'karyotype.txt' from circos_input_karyotype_chr.first()
        file 'genes.txt' from circos_input_genes_chr.first()
        file 'gaps.txt' from circos_input_gaps_chr.first()
        file 'core.txt' from core_comp_circos_chr.first()
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

    bin_target_mapping.splitCsv(sep: "\t").set { circos_binmap }
    circos_binconffile = file(params.CIRCOS_BIN_CONFIG_FILE)

    process circos_run_bin {
        // this process can fail
        errorStrategy 'ignore'

        input:
        file 'links.txt' from circos_input_links_bin.first()
        file 'karyotype.txt' from circos_input_karyotype_bin.first()
        file 'genes.txt' from circos_input_genes_bin.first()
        file 'gaps.txt' from circos_input_gaps_bin.first()
        file 'core.txt' from core_comp_circos_bin.first()
        val circos_binconffile
        val binmap from circos_binmap

        output:
        set file('bin.png'), val(binmap) into circos_bin_output

        """
        circos  -conf ${circos_binconffile} -param image/file=bin.png  \
                -param chromosomes='${binmap[0]};${binmap[1]}' \
                -param chromosomes_reverse=${binmap[0]} \
                -param chromosomes_scale='${binmap[0]}:0.4r'
        """
    }

    circos_output.subscribe {
        if (params.print_paths) {
          println it[0]
        }
        if (params.dist_dir) {
          it[0].copyTo(params.dist_dir + "/chr" + it[1][0] + ".png")
        }
    }

    circos_bin_output.subscribe {
        if (params.print_paths) {
          println it[0]
        }
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
        file 'pseudo.fasta.gz' into embl_full_seq

        """
        gt inlineseq_add -seqfile pseudo.fasta.gz -matchdescstart -force \
          -o full.gff3 pseudo.gff3
        """
    }

    process make_embl {
        afterScript 'rm -rf 1 1.*'

        input:
        file 'embl_in.gff3' from embl_full_gff
        file embl_full_seq
        val go_obo

        output:
        file '*.embl' into embl_out

        script:
        if (params.embl_ena_submission)
            """
            zcat ${embl_full_seq} > 1 && gff3_to_embl.lua -e -o embl_in.gff3 ${go_obo} '${params.EMBL_ORGANISM}' 1
            """
        else
            """
            zcat ${embl_full_seq} > 1 && gff3_to_embl.lua -o embl_in.gff3 ${go_obo} '${params.EMBL_ORGANISM}' 1
            """
    }

    embl_out.flatMap().subscribe {
        if (params.print_paths) {
          println it
        }
        if (params.dist_dir) {
          it.copyTo(params.dist_dir)
        }
    }
}

// REPORT CREATION
// ===============

specfile = file(params.SPECFILE)
process make_report {
    validExitStatus 0,1,2

    input:
    set file('pseudo.fasta.gz'), file('scaf.fasta.gz') from report_inseq
    set file('pseudo.gff3'), file('scaf.gff3') from report_gff3
    val specfile

    output:
    set file('pseudo.report.html'), file('scaf.report.html') into report_output

    """
    touch pseudo.report.html
    touch scaf.report.html
    gt speck -specfile ${specfile} -matchdescstart -seqfile pseudo.fasta.gz \
      -provideindex -typecheck so -output "${params.SPECK_TEMPLATE}" \
      pseudo.gff3 > pseudo.report.html || true
    gt speck -specfile ${specfile} -matchdescstart -seqfile scaf.fasta.gz \
      -provideindex -typecheck so -output "${params.SPECK_TEMPLATE}" \
      scaf.gff3 > scaf.report.html || true
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
    if (params.print_paths) {
      println it
    }
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
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}

out_seq.subscribe {
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}

result_agp.subscribe {
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}

result_gaf.collectFile().subscribe {
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

result_ortho.collectFile().subscribe {
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

result_protein.subscribe {
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

stats_output.subscribe {
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      it.copyTo(params.dist_dir)
    }
}

report_output.subscribe {
    if (params.print_paths) {
      println it
    }
    if (params.dist_dir) {
      for (file in it) {
        file.copyTo(params.dist_dir)
      }
    }
}
