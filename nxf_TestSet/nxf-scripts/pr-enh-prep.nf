nextflow.preview.dsl=2


ch_enhancer_inp = Channel.fromPath( params.genomes[params.refgen].enhancer_bed )
ch_promoter_inp = Channel.fromPath( params.genomes[params.refgen].promoter_gtf )

if (params.refgen=="hg38" || params.refgen=="mm10"){
  promoter_processing = "pt.process_promoter_v2"
}else if(params.refgen=="hg19"){
  promoter_processing = "pt.process_promoter_v1"
}else{
    exit 1, "supported ref genome hg19, hg38, mm10"
}

enhBed_headers=params.genomes[params.refgen].enhBed_headers
enh_window = params.enhancer_window.toInteger() + params.augment_length.toInteger()
pr_window = params.promoter_window.toInteger() + params.augment_length.toInteger()

// promoter pre-processing
process PREPROCESS_PROMOTER {
    publishDir "${params.publish_dir}"

    input:
    path(input)

    output:
    // stdout result
    tuple file("promoter.${params.refgen}.bed"), file("promoter_bg.${params.refgen}.bed")

    script:
    """
    python -c "import preptools as pt; \
    $promoter_processing('$input', \
        $pr_window, \
        $params.all_chrom, \
        'promoter.${params.refgen}.bed')"

    python -c "import preptools as pt; \
    $promoter_processing('$input', \
        $params.bgWindow, \
        $params.all_chrom, \
        'promoter_bg.${params.refgen}.bed')"
    """
    // func( allfield_bed, bg_path, headers, window=bg_window)
}

// enhancer pre-processing
process PREPROCESS_ENHANCER {
    publishDir "${params.publish_dir}"

    input:
    path(input)

    output:
    // stdout result
    tuple file("enhancer.${params.refgen}.bed"), file("enhancer_bg.${params.refgen}.bed")


    script:
    """
    python -c "from preptools import process_enhancer; \
    process_enhancer('$input', \
        $enh_window, \
        $params.all_chrom, \
        $enhBed_headers, \
        'enhancer.${params.refgen}.bed')"

    python -c "from preptools import process_enhancer; \
    process_enhancer('$input', \
        $params.bgWindow, \
        $params.all_chrom, \
        $enhBed_headers, \
        'enhancer_bg.${params.refgen}.bed')"
    """
}

process BRUTE_FORCE_PR_ENH_PAIRS {
    publishDir "${params.publish_dir}"

    input:
    tuple path(enh_bed), path(pr_bed)

    output:
    path("pairs_bruteForce.${params.refgen}.csv")

    script:
    """
    #!/usr/bin/env bash
    gen_enhPr_interac.sh $params.refgen $params.max_dist_pr_enh $enh_window $pr_window $enh_bed $pr_bed
    """
}

workflow {
    println "${params.publish_dir}"
    ch_enh_out = PREPROCESS_ENHANCER( ch_enhancer_inp )
    ch_pr_out = PREPROCESS_PROMOTER( ch_promoter_inp )

    ch_enh_pr = ch_enh_out.map{ it -> it[0]}.combine(ch_pr_out.map{ it -> it[0]})
    ch_bruteforce_enh_pr = BRUTE_FORCE_PR_ENH_PAIRS(ch_enh_pr)

    ch_enh_out.view()
    ch_pr_out.view()
    ch_bruteforce_enh_pr.view()

}