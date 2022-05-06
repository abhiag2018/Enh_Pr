nextflow.preview.dsl=2

ch_input_fasta = Channel.fromPath(params.genomes[params.refgen].fasta)
ch_input_fasta_fai = Channel.fromPath(params.genomes[params.refgen].fasta_idx)


ch_regelem = Channel.fromPath("$params.regelem_inp")
    .splitCsv(header:true, sep:',')
    .map { row -> ["$row.regelem", file("$params.store_dir/$row.main", checkIfExists: true) ] }
    .filter{ it[0] == "enhancer" || it[0] == "promoter" }
    .map{ it -> it[1] }


// if (params.dev){
//     ch_enhancer_bed_prep = Channel.fromPath("enhancer_mm.dev.bed")
//     ch_promoter_bed_prep = Channel.fromPath("promoter_mm.dev.bed")
// }else{
//     Channel.fromPath("enh_bed.csv")
//       .splitCsv(header:true, sep:',')
//       .map { row -> file("$params.store_dir/$row.enhancer", checkIfExists: true)  }
//       .set {ch_enhancer_bed_prep}

//     Channel.fromPath("pr_bed.csv")
//       .splitCsv(header:true, sep:',')
//       .map { row -> file("$params.store_dir/$row.promoter", checkIfExists: true) }
//       .set {ch_promoter_bed_prep}    
// }

process GEN_FASTA{
    input:
    tuple path(reg_elem_bed), path(fasta), path(fasta_fai)

    output:
    file("${reg_elem_bed.baseName}.fa")

    script:
    """
    bedtools getfasta -fi $fasta -bed $reg_elem_bed -name -fo ${reg_elem_bed.baseName}.fa
    """
}


process GEN_CSV{
    input:
    path(reg_elem_fasta)

    output:
    file("${reg_elem_fasta.baseName}.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("Biostrings")

    fastaFile <- readDNAStringSet("$reg_elem_fasta")
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    df <- data.frame(seq_name, sequence)

    write.csv(df,"${reg_elem_fasta.baseName}.csv", row.names = FALSE)
    """

}


process GEN_NPZ{
    publishDir "${params.publish_dir}"

    input:
    path(reg_elem_fa_csv)

    output:
    file("${reg_elem_fa_csv.baseName}.dnaseq.npz")

    script:
    """
    #!/usr/bin/env python
    from sklearn.preprocessing import OneHotEncoder
    import pandas as pd
    import numpy as np
    import csv
    import h5py
    enc = OneHotEncoder(handle_unknown='ignore')
    init = np.array(list('ATCG')).reshape(4,1)
    enc.fit(init)

    regElm_fa_df = pd.read_csv("${reg_elem_fa_csv}")

    regElm_onehot = regElm_fa_df.sequence.apply(lambda line:enc.transform(np.array(list(line.upper())).reshape(len(line),1)).toarray())
    regElm_fa_df.sequence = regElm_onehot.apply(lambda x:list(np.array(x,dtype=int).reshape(x.shape[0]*x.shape[1])))
    X = regElm_fa_df.apply(lambda df:df.seq_name.split("::"), axis=1,result_type ='expand')
    regElm_fa_df['name'] = X[0]
    regElm_fa_df['loc'] = X[1]
    # regElm_fa_df[['sequence','name','loc']].to_csv("${reg_elem_fa_csv.baseName}.dnaseq.csv",index=False)
    # with h5py.File("${reg_elem_fa_csv.baseName}.dnaseq.h5", 'w') as hf:
    #     hf.create_dataset('sequence',  data=np.stack(regElm_fa_df['sequence']))
    #     hf.create_dataset('name',  data=regElm_fa_df['name'], dtype=h5py.special_dtype(vlen=str))
    #     hf.create_dataset('loc',  data=regElm_fa_df['loc'], dtype=h5py.special_dtype(vlen=str))
    # # out = h5py.File("h5-file",'r')['sequence']
    np.savez("${reg_elem_fa_csv.baseName}.dnaseq.npz",sequence=np.stack(regElm_fa_df['sequence']),name=regElm_fa_df['name'],loc=regElm_fa_df['loc'])
    """
}


// workflow genseq_promoter{ emit: GEN_FASTA(ch_promoter_bed_prep,ch_input_fasta,ch_input_fasta_fai) | GEN_CSV | GEN_NPZ }
// workflow genseq_enhancer{ emit: GEN_FASTA(ch_enhancer_bed_prep,ch_input_fasta,ch_input_fasta_fai) | GEN_CSV | GEN_NPZ }

workflow{
  ch_input = ch_regelem.combine(ch_input_fasta).combine(ch_input_fasta_fai)
  ch_output = GEN_FASTA(ch_input) | GEN_CSV | GEN_NPZ 
  ch_output.view()

  // genseq_promoter().view()
  // genseq_enhancer().view()
}












