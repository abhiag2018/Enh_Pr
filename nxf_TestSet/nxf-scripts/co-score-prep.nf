nextflow.preview.dsl=2
// include { promoter_bed; enhancer_bed } from "${params.codeDir}/NN-feature-prep/modules/pr-enh-prep"
include {CHROM_OPENN_SCORE_PROFILE as PROFILE_ENH} from './module'
include {CHROM_OPENN_SCORE_PROFILE as PROFILE_PR} from './module'

ch_chrom_val = Channel.fromList(Eval.me(params.species[params.refgen.substring(0,2)].all_chrom))
ch_fasta_idx = Channel.value( params.genomes[params.refgen].fasta_idx )

Channel.fromPath("$params.regelem_inp")
    .splitCsv(header:true, sep:',')
    .filter{ "$it.regelem" == "enhancer" || "$it.regelem" == "promoter" }
    .map { row -> ["$row.regelem", "$params.store_dir/$row.main", "$params.store_dir/$row.background" ] }
    .map{ it -> [it[0], file(it[1], checkIfExists: true), file(it[2], checkIfExists: true)] }
    .set{ ch_regelem }

ch_regelem
  .filter{ it[0] == "enhancer" }
  .map{ it -> [it[1], it[2]] }
  .set{ ch_enhancer_bed_prep }

ch_regelem
  .filter{ it[0] == "promoter" }
  .map{ it -> [it[1], it[2]] }
  .set{ ch_promoter_bed_prep }

ch_input_bam = Channel.fromPath("$params.bamInput", checkIfExists: true)
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, file(row.bam, checkIfExists: true) ]  }
if (params.dev){
  ch_input_bam = ch_input_bam.filter{ it[0] =~ ~/type66/ }
}


// DNase/ATAC preprocessing step 1
// generate .bam.bai index file from DNase/ATAC-seq .bam files
process GEN_BAM_INDEX {
    input:
    tuple val(cell), val(rep), path(bamfile)

    output:
    tuple val(cell), val(rep), file("${bamfile}.bai"), file("${bamfile}")

    script:  
    """
      samtools index $bamfile ${bamfile}.bai
    """
}

// DNase/ATAC preprocessing step 2
// chromatin openness score construction from .bam file:
// at each genomic location count the number of reads intersecting with that location
// .. for each .bam and for each chromosome
process CHROM_OPENN_SCORE {
    input:
    tuple val(cell), val(rep), path(bamfile_index), path(bamfile), val(chrom)
    val(fasta_idx)

    output:
    tuple val(cell), val(rep), path("*.npz")


    script:
    """
    #!/usr/bin/env python
    from preptools import getBamCounts
    import pandas as pd
    df = pd.read_csv("$fasta_idx",sep="\t",header=None)
    len_val=dict(zip(df[0], df[1]))["${chrom}"]
    getBamCounts('${bamfile}', 
        '${chrom}', 
        len_val, 
        outputf='${cell}.rep${rep}.${chrom}.npz')
    """
}

workflow {
    ch_bam_index = GEN_BAM_INDEX(ch_input_bam) 
    ch_bam_chrom = ch_bam_index.combine(ch_chrom_val)
    ch_bam_chrom_readcounts = CHROM_OPENN_SCORE(ch_bam_chrom,ch_fasta_idx)

    ch_bam_readcounts_grouped = ch_bam_chrom_readcounts.groupTuple(by: [0,1])

    PROFILE_ENH(ch_bam_readcounts_grouped.combine(ch_enhancer_bed_prep).combine(Channel.value("enhancer")))
        .mix(PROFILE_PR(ch_bam_readcounts_grouped.combine(ch_promoter_bed_prep).combine(Channel.value("promoter"))))
        .view()
}











