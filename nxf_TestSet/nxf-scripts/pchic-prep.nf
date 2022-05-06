// with CONNECTION
import groovy.json.JsonSlurper
nextflow.preview.dsl=2
include {SPLIT_BY_CHROM; SPLIT_HIC; MATCH_HIC_REGELEM; SELECT_UNIQ_REGELEM; COMBINE_HIC; INTEGRATE; AUGMENT_DATA; COMBINE_FILES} from './module'

include {KNAPSACK as KNAPSACK_0} from './module'

include {SUBLIST as SUBLIST_0} from './module'
include {SUBLIST as SUBLIST_1} from './module'
include {SUBLIST as SUBLIST_2} from './module'
include {SUBLIST as SUBLIST_3} from './module'

ch_cellTypes = Channel.fromList(params.cellTypes)

Channel.fromPath("$params.regelem_inp")
    .splitCsv(header:true, sep:',')
    .map { row -> ["$row.regelem", "$params.store_dir/$row.main" ] }
    .filter{ it[0] == "brutePairs" }
    .map{ it -> file(it[1], checkIfExists: true) }
    .set {ch_hic_all}


Channel.fromPath("$params.dnaseq_inp")
  .splitCsv(header:true, sep:',')
  .map { row -> [ file("$params.store_dir/$row.enhancer", checkIfExists: true), file("$params.store_dir/$row.promoter", checkIfExists: true) ]  }
  .set {ch_gen_seq}

workflow prep_pchic{ 
    // cheeky_hugle //all
    // determined_mestorf // split = [0.2, 0.2] // test / val / train

    ch_hic_input = Channel.fromPath(params.hic_input)
        .splitCsv(header:true, sep:',')
        .map { row -> ["$row.hic_input_name", "$row.path" ] }
        .map{ it -> [it[0], file(it[1], checkIfExists: true)] }
    ch_hic_input_list = ch_hic_input.map{ it -> it[0] }

    ch_hic_input_chrom_0 = SPLIT_BY_CHROM(ch_hic_input).transpose()
    ch_hic_input_chrom = ch_hic_input_chrom_0.map{ it -> [it[0], (it[1] =~ "(chr.*?)\\.")[0][1], it[1]] }
    if (params.dev){
        ch_hic_input_chrom = ch_hic_input_chrom
            .splitText(by:1000,keepHeader:true,file:true)
            .unique{ [it[0],it[1]] }
    }

    ch_hic_cell_regElements_0 = SPLIT_HIC(ch_hic_input_chrom).transpose().combine(ch_cellTypes)

    ch_hic_cell_regElements = ch_hic_cell_regElements_0.filter{ file(it[2]).countLines()>1 }

    ch_hic_match = MATCH_HIC_REGELEM(ch_hic_input_list.combine(ch_gen_seq).combine(ch_hic_cell_regElements, by: 0)) | SELECT_UNIQ_REGELEM
    ch_hic_match_uniq = ch_hic_match.groupTuple(by: [0,1,2]) | COMBINE_HIC
    _ch_hic_match_uniq_all = INTEGRATE(ch_hic_match_uniq.combine(ch_hic_all)) | AUGMENT_DATA


    // separate data *by chromosome* for cross-validation
    // we have to partition the chromosomes into #params.split parts 
    // I use knapsack algorithm to partition the chromosomes in #params.split subsets 


    ch_list_frac = params.split+[1-params.split.sum()]


    _ch_hic_match_uniq_all
        // .filter{ it[3] =~  partition_regex[i] }
        .map{it -> [it[0],it[1],it[2],it[3],file(it[3]).countLines(),it[4],it[5]]} // count lines
        .groupTuple(by: [0,2]) // combine chromosome data for the same hic-file and cell type
        .map{ it -> [[it[0],it[1],it[2],it[3],it[5],it[6]],it[4], ch_list_frac]}//.cross(ch_list_frac) 
        .set{ __ch_hic_match_uniq_all }
    KNAPSACK_0(__ch_hic_match_uniq_all) // knapsack to find the best cross-validation partition
        .map{ it -> [it[0], file(it[1], checkIfExists: true)]}
        // .map{ it -> [it[0][0],it[0][1],it[0][2],it[0][3],it[1],it[0][4],it[0][5]]}
        // .map{ it -> [ [it[0],it[1],it[2],it[3],it[5],it[6]], file(it[4], checkIfExists: true)] }
        .splitCsv(elem:1, header:true, sep:';')
        .map { it -> [it[0], it[1]['num_cols'], it[1]['label'], it[1]['list_index'] ] }
        .set{ ___ch_hic_match_uniq_all }    


    ___ch_hic_match_uniq_all
        .map{ it -> [[it[0][0],it[0][2],it[2],it[1]], it[0][1],it[3]] }
        .set{ ch_chr }


    ___ch_hic_match_uniq_all
        .map{ it -> [[it[0][0],it[0][2],it[2],it[1]], it[0][3],it[3]] }
        .set{ ch_chr_intgr }

    ___ch_hic_match_uniq_all
        .map{ it -> [[it[0][0],it[0][2]],it[0][4],"${it[0][0]}.${it[0][2]}.pos.csv"] }
        .unique()
        .set{ ch_pos }

    ___ch_hic_match_uniq_all
        .map{ it -> [[it[0][0],it[0][2]],it[0][5],"${it[0][0]}.${it[0][2]}.unmatched.csv"] }
        .unique()
        .set{ ch_unmat }

    ch_hic_match_uniq_aug = SUBLIST_0(ch_chr) // select chromosome subset computed from KNAPSACK

    ch_hic_match_uniq_aug_1 = SUBLIST_1(ch_chr_intgr) // select chromosome subset computed from KNAPSACK
    ch_hic_match_uniq_aug = ch_hic_match_uniq_aug.join(ch_hic_match_uniq_aug_1, by:0, remainder:true)
        .map{it -> [[it[0][0],it[0][1],it[0][2],it[0][3],it[1]],it[2],"${it[0][0]}.${it[0][1]}.${it[0][2]}.intgr.csv"] }

    ch_hic_match_uniq_aug.collectFile(storeDir: "${params.publish_dir}") {it -> ["${it[0][0]}.${it[0][1]}.${it[0][2]}.README.txt", "lines: ${it[0][3]}"+'\n'+"chrom: ${it[0][4]}"+'\n'] }

    ch_hic_match_uniq_aug_combined = ch_hic_match_uniq_aug.mix(ch_pos,ch_unmat) | COMBINE_FILES

    emit:
    // ch_hic_match_uniq_aug
    ch_hic_match_uniq_aug_combined

}

workflow {
    ch_out = prep_pchic()
    ch_out.view()
    ch_out.count().view()
}








