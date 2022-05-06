process CONSOLIDATE_OUT {
    stageInMode 'link'
    publishDir "$params.publish_dir"
    tag "data.${cellType}.h5"
    label 'process_intense'

    input:
    tuple val(cellType), val(label), path(data_h5_list), val(app_str)

    output:
    tuple val(cellType), val(label), path("data.${cellType}.rep${label}${app_str}.h5")

    script:
    """
    #!/usr/bin/env python

    import h5py
    import pandas as pd
    import numpy as np
    from natsort import natsorted

    data_h5_list = natsorted("$data_h5_list".split())
    dset_size = 0
    for h5 in data_h5_list:
        with h5py.File(h5,'r') as f:
            dset_size += f['label'].shape[0]

    with h5py.File("data.${cellType}.rep${label}${app_str}.h5",'w') as h5f0:
        with h5py.File(data_h5_list[0],'r') as h5f:
            h5f0_dset = {}
            for k in h5f.keys():
                h5f0_dset[k] = h5f0.create_dataset(k, (dset_size,)+h5f[k].shape[1:], dtype = h5f[k].dtype)

        ind=0
        for h5 in data_h5_list:
            with h5py.File(h5,'r') as h5f:
                len_chunk = h5f['label'].shape[0]

                for k in h5f0.keys():
                    h5f0_dset[k].write_direct(h5f[k][:], dest_sel=np.s_[ind:ind+len_chunk])

                ind += len_chunk
    """
}

// pchic-prep.nf
all_chrom = params.species[params.refgen.substring(0,2)].all_chrom
process SPLIT_BY_CHROM{
    memory '10 GB'

    input:
    tuple val(hic_input_name), path(input_csv)

    output:
    tuple val(hic_input_name), path("${input_csv.baseName}.chr*.tsv")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    pchicDF = pd.read_csv("$input_csv", sep="\t",dtype={'baitChr':str,'oeChr':str})
    for chrom in $all_chrom:
        chrom = chrom[3:]
        pchicDF[(pchicDF.oeChr.apply(str)==chrom) & (pchicDF.baitChr.apply(str)==chrom)].to_csv("${input_csv.baseName}.chr"+chrom+".tsv",header=True,sep="\t",index=False)
    pchicDF[pchicDF.oeChr != pchicDF.baitChr].to_csv("${input_csv.baseName}.chr"+"CROSS"+".tsv",header=True,sep="\t",index=False)
    """
}

hic_split_num = params.dev ? 2 : params.hic_split_process

// PCHi-C processing : step 1
// split input tsv file for PCHi-C data into params.hic_split_process parts
process SPLIT_HIC {
    input:
    tuple val(hic_input_name), val(chrom), path(hic_input)

    output:
    tuple val(hic_input_name), val(chrom), path("${hic_input.baseName}.*.tsv")

    script:
    """
    python -c "import preptools as pt; \
    pt.splitCSV('$hic_input', \
        [], \
        readArgs = {'delimiter':'\t','dtype':{'baitChr':str,'oeChr':str}}, \
        writeArgs = {'header':True,'sep':'\t','index':False}, prefix='${hic_input.baseName}.', \
        split_num=$hic_split_num, \
        suffix='.tsv')"
    """

}

min_threshold = 0 //threshold for selecting all values 
//params.pos_threshold

// PCHi-C processing : step 2
// match the elements in the tsv file to elements in the regulatory elements lists
// 6 min
process MATCH_HIC_REGELEM {
    memory { 30.GB }

    input:
    tuple val(hic_input_name), path(enhancer_dnaseq_gz), path(promoter_dnaseq_gz), val(chrom), path(hic_input), val(cellType)

    output:
    tuple val(hic_input_name), val(chrom), val(cellType), path('*.pkl')

    script:  
    ext = enhancer_dnaseq_gz.extension
    if (ext=="gz") {
        enhancer_npz = enhancer_dnaseq_gz.baseName
        promoter_npz = promoter_dnaseq_gz.baseName
    }
    else{
        enhancer_npz = enhancer_dnaseq_gz
        promoter_npz = promoter_dnaseq_gz
    }
    """
    if [ "$ext" == ".gz" ]; then
        gzip -df $enhancer_dnaseq_gz
        gzip -df $promoter_dnaseq_gz
    fi
    python -c "import preptools as pt; \
    import os; \
    pt.concat_PCHiC_PE('$hic_input', \
        '$promoter_npz', \
        '$enhancer_npz', \
        selectCell='$cellType', \
        threshold = ${min_threshold}, \
        outputF='.'.join('$hic_input'.split('.')[:-1])+'.matched_${cellType}.pkl', \
        sampleFrac=None)
    os.system('rm $promoter_npz')
    os.system('rm $enhancer_npz')"
    """ 
}

enh_window = params.enhancer_window.toInteger() + params.augment_length.toInteger()
pr_window = params.promoter_window.toInteger() + params.augment_length.toInteger()

// PCHi-C processing : step 4
// select unique promoter/enhancer from matched elements
process SELECT_UNIQ_REGELEM {
    memory '100 GB'
    time '3d'
    errorStrategy 'finish'

    input:
    tuple val(hic_input_name), val(chrom), val(cellType), path(matched_hic_parts)

    output:
    tuple val(hic_input_name), val(chrom), val(cellType), path("${matched_hic_parts.baseName}.uniq.csv")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np

    cols = ['baitChr', 'baitStart', 'baitEnd', 'oeChr', 'oeStart', 'oeEnd',"$cellType"]
    cols_1=['enhancer_chrom', 'enhancer_start', 'enhancer_end', 'enhancer_name', 'promoter_chrom', 'promoter_start', 'promoter_end', 'promoter_name','label',"$cellType"]
    pchicDF = pd.read_pickle("${matched_hic_parts}")
    # pchicDF = pd.concat([pd.read_pickle(inp) for inp in "${matched_hic_parts}".split()])

    baitprF = pchicDF['baitPr'].apply(lambda x:np.array(x).size>0)
    baitenhF = pchicDF['baitEnh'].apply(lambda x:np.array(x).size>0)
    oeprF = pchicDF['oePr'].apply(lambda x:np.array(x).size>0)
    oeenhF = pchicDF['oeEnh'].apply(lambda x:np.array(x).size>0)

    pehic = pchicDF[baitprF&oeenhF][cols+['baitPr','oeEnh']].copy()
    ephic = pchicDF[baitenhF&oeprF][cols+['oePr','baitEnh']].copy()

    ## select those interactions where the bait intersects with a unique gene name
    unique_gene_pe = pehic.baitPr.apply(lambda X:np.unique(list(map(lambda x:x.split(";")[1],np.array(X)[:,1])))).apply(lambda x:len(x)==1)
    pehic = pehic[unique_gene_pe].copy()

    if len(pehic)>0:
        ## now, out of the unique gene name bait, oe interactions select one promoter/enhancer
        pehic['uniqPr'] = pehic.apply(lambda rw:np.array(list(map(int,np.array(rw['baitPr'])[:,0])))-(rw.baitStart+rw.baitEnd)//2, axis=1).apply(abs).apply(np.argmin)
        pehic['baitPr'] = pehic.apply(lambda rw:np.array(rw['baitPr'])[rw.uniqPr,:],axis=1)
        pehic['uniqEnh'] = pehic.apply(lambda rw:np.array(list(map(int,np.array(rw['oeEnh'])[:,0])))-(rw.oeStart+rw.oeEnd)//2, axis=1).apply(abs).apply(np.argmin)
        pehic['oeEnh'] = pehic.apply(lambda rw:np.array(rw['oeEnh'])[rw.uniqEnh,:],axis=1)

        ## rename/prettify & select useful columns 
        pehic['promoter_chrom'] = 'chr'+pehic['baitChr']
        pehic['promoter_start'] = pehic['baitPr'].apply(lambda x:int(x[0])-$pr_window//2)
        pehic['promoter_end'] = pehic['baitPr'].apply(lambda x:int(x[0])+$pr_window//2)
        pehic['promoter_name'] = pehic['baitPr'].apply(lambda x:x[1])
        pehic['enhancer_chrom'] = 'chr'+pehic['oeChr']
        pehic['enhancer_start'] = pehic['oeEnh'].apply(lambda x:int(x[0])-$enh_window//2)
        pehic['enhancer_end'] = pehic['oeEnh'].apply(lambda x:int(x[0])+$enh_window//2)
        pehic['enhancer_name'] = pehic['oeEnh'].apply(lambda x:x[1])
        pehic['label'] = 1
        pehic = pehic[cols_1]
    else:
        pehic=pd.DataFrame(columns=cols_1)

    ## select those interactions where the oe intersects with a unique gene name
    unique_gene_ep = ephic.oePr.apply(lambda X:np.unique(list(map(lambda x:x.split(";")[1],np.array(X)[:,1])))).apply(lambda x:len(x)==1)
    ephic = ephic[unique_gene_ep].copy()

    if len(ephic)>0:
        ephic['uniqPr'] = ephic.apply(lambda rw:np.array(list(map(int,np.array(rw['oePr'])[:,0])))-(rw.oeStart+rw.oeEnd)//2, axis=1).apply(abs).apply(np.argmin)
        ephic['oePr'] = ephic.apply(lambda rw:np.array(rw['oePr'])[rw.uniqPr,:],axis=1)
        ephic['uniqEnh'] = ephic.apply(lambda rw:np.array(list(map(int,np.array(rw['baitEnh'])[:,0])))-(rw.baitStart+rw.baitEnd)//2, axis=1).apply(abs).apply(np.argmin)
        ephic['baitEnh'] = ephic.apply(lambda rw:np.array(rw['baitEnh'])[rw.uniqEnh,:],axis=1)


        ephic['promoter_chrom'] = 'chr'+ephic['oeChr']
        ephic['promoter_start'] = ephic['oePr'].apply(lambda x:int(x[0])-$pr_window//2)
        ephic['promoter_end'] = ephic['oePr'].apply(lambda x:int(x[0])+$pr_window//2)
        ephic['promoter_name'] = ephic['oePr'].apply(lambda x:x[1])
        ephic['enhancer_chrom'] = 'chr'+ephic['baitChr']
        ephic['enhancer_start'] = ephic['baitEnh'].apply(lambda x:int(x[0])-$enh_window//2)
        ephic['enhancer_end'] = ephic['baitEnh'].apply(lambda x:int(x[0])+$enh_window//2)
        ephic['enhancer_name'] = ephic['baitEnh'].apply(lambda x:x[1])
        ephic['label'] = 1
        ephic = ephic[cols_1]
    else:
        ephic=pd.DataFrame(columns=cols_1)

    print("duplicated rows:")
    print("\tpehic: %s/%s"%(np.sum(pehic.duplicated()),len(pehic)))
    print("\tephic: %s/%s"%(np.sum(ephic.duplicated()),len(ephic)) )

    pd.concat([ephic,pehic]).to_csv("${matched_hic_parts.baseName}.uniq.csv", index=False)
    # pd.concat([ephic,pehic]).to_csv("${chrom}.${cellType}.uniq.csv", index=False)
    """
}

// PCHi-C processing : step 3
// combine matched hic .tsv files and drop duplicates
process COMBINE_HIC {
    input:
    tuple val(hic_input_name), val(chrom), val(cellType), path(chrom_hic_csv)

    output:
    tuple val(hic_input_name), val(chrom), val(cellType), path("*.csv")

    script:
    if (chrom instanceof List){
        chrom="all"
    }
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    m_files = "${chrom_hic_csv}".split()
    out_fname = '.'.join(m_files[0].split('.')[:-5])+".${cellType}.${chrom}.csv"
    df_output = pd.concat([pd.read_csv(inp) for inp in m_files])

    
    cols = df_output.columns[:-1]
    df_output_ = df_output.groupby(list(cols)).apply(lambda df: df[df['$cellType'] == df['$cellType'].max()] ).reset_index(drop=True).drop_duplicates(ignore_index=True)

    print("duplicated rows:")
    print( "\tdf_output: %s/%s"%( df_output.shape[0] - df_output_.shape[0], len(df_output) ) )
    # df_output = df_output.drop_duplicates(ignore_index=True)

    df_output_.to_csv(out_fname, index=False)
    """
}

// integrate the filtered positive interactions with the 
// list of all possible pr-enh interactions for the species
process INTEGRATE{
    memory '30 GB'
    time '3d'
    errorStrategy 'finish'

    input:
    tuple val(hic_input_name), val(chrom), val(cellType), path(hic_cell), path(hic_all)

    output:
    tuple val(hic_input_name), val(chrom), val(cellType), path("${hic_cell.baseName}.intgr.csv"), path(hic_cell), path("${hic_input_name}.unmatched.${cellType}.${chrom}.txt")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd, numpy as np

    hic_all = pd.read_csv("$hic_all")
    hic_chr = hic_all[hic_all["enhancer_chrom"]=="$chrom"]
    pat_df = pd.read_csv("$hic_cell")

    inter = hic_chr.merge(pat_df, how='left', on=list(hic_chr.columns), indicator=True)
    inter.label = inter._merge.apply(lambda x: 1 if x=="both" else 0)
    inter[pat_df.columns].to_csv("${hic_cell.baseName}.intgr.csv",index=False)


    unmat = pat_df.merge(hic_chr, how='left', on=list(hic_chr.columns), indicator=True)
    unmat[unmat._merge == "left_only"][list(pat_df.columns)].to_csv("${hic_input_name}.unmatched.${cellType}.${chrom}.txt",index=False)
    """
}

// augment datasets if needed for generating training set
// also set the enhancer window (enh_window) and promoter window (pr_window) correctly
process AUGMENT_DATA{
    memory '30 GB'
    time '3d'
    errorStrategy 'finish'

    input:
    tuple val(hic_input_name), val(chrom), val(cellType), path(hic_out), path(hic_cell), path(unmatched_rows)

    output:
    tuple val(hic_input_name), val(chrom), val(cellType), path("${hic_out.baseName}.*.csv"), path(hic_cell), path(unmatched_rows) mode flatten

    script:
    """
    #!/usr/bin/env python
    import pandas as pd, numpy as np; np.random.seed(42)
    import random; random.seed(21)
    import math
    import string

    enh_steps = np.arange(0,$params.augment_length,$params.augment_step)
    pr_steps = np.arange(0,$params.augment_length,$params.augment_step)

    enh_window = $params.enhancer_window
    pr_window = $params.promoter_window

    neg_to_pos = math.inf if "$params.neg_pos_interac_ratio" == 'inf' else $params.neg_pos_interac_ratio

    pos_threshold = $params.pos_threshold

    df0 = pd.read_csv("$hic_out")

    split = []
    # split = $params.split
    # test_val_train = ['test','val','train']
    test_val_train = list(string.ascii_uppercase)[:(len(split)+1)]
    labels = ['all'] if len(split)==0 else [x+str(y).replace('.','_') for x,y in zip(test_val_train,split + [1.0-sum(split)])]

    pos_df = df0[df0.label==1].sample(frac=1.0)
    neg_df = df0[df0.label==0].sample(frac=1.0)
    if len(split)>0:
        N_pos = [int(x*pos_df.shape[0]) for x in split]
        N_pos = N_pos + [pos_df.shape[0]-sum(N_pos)]
        N_neg = [int(x*neg_df.shape[0]) for x in split]
        N_neg = N_neg + [neg_df.shape[0]-sum(N_neg)]
        
        df_list_pos = []
        prev_rows=0
        for r in N_pos:
            df_list_pos += [pos_df.iloc[prev_rows:prev_rows+r].copy(deep=True)]
            prev_rows += r

        df_list_neg = []
        prev_rows=0
        for r in N_neg:
            df_list_neg += [neg_df.iloc[prev_rows:prev_rows+r].copy(deep=True)]
            prev_rows += r

        df_list = []
        for pos_, neg_ in zip(df_list_pos,df_list_neg):
            df_list += [pd.concat([pos_, neg_])]
    else:
        df_list = [df0]
    
    def aug_random(df, size_multip):
        if size_multip==0:
            _df = df[0:0]
        elif size_multip == int(size_multip):
            _df = pd.concat([df]*int(size_multip), ignore_index=True)
        elif size_multip>1:
            frac_mult = size_multip-int(size_multip)
            _df_1 = pd.concat( [df]*int(size_multip), ignore_index=True)
            _df_2 = df.sample(n=int(len(df)*frac_mult))
            _df = pd.concat([_df_1, _df_2], ignore_index=True)
        else:
            _df = df.sample(n=int(len(df)*size_multip))

        if size_multip>1:
            _df.enhancer_start = _df.enhancer_start.apply(lambda x:x+random.choice(enh_steps) )
            _df.promoter_start = _df.promoter_start.apply(lambda x:x+random.choice(pr_steps) )
        else:
            _df.enhancer_start = (_df.enhancer_start+ _df.enhancer_end)//2 - enh_window//2
            _df.promoter_start = (_df.promoter_start+ _df.promoter_end)//2 - pr_window//2

        _df.enhancer_end = _df.enhancer_start + enh_window
        _df.promoter_end = _df.promoter_start + pr_window

        return _df


    for i, df0 in enumerate(df_list):
        pos_mult_fac = $params.hic_augment_factor if i == len(df_list)-1 else 1

        pos_df = df0[df0.label==1]
        pos_df = aug_random(pos_df[pos_df["$cellType"]>pos_threshold], pos_mult_fac)

        neg_df = df0[df0.label==0]
        if len(neg_df) == 0 & len(pos_df) != 0:
            raise NameError("no neg interactions!!")
        if neg_to_pos != math.inf:
            neg_mult = 1 if len(neg_df)==0 else len(pos_df)/len(neg_df)*neg_to_pos
            neg_df = aug_random(neg_df, neg_mult)
        pd.concat( [pos_df, neg_df], ignore_index=True).to_csv(f"${hic_out.baseName}.{labels[i]}.csv", index=False)
    # .to_csv("${hic_out.baseName}.train.csv", index=False)
    """
}


// process COMBINE_CHR{
//     publishDir "${params.publish_dir}"

//     input:
//     tuple val(identifier), val(hic_input_name), val(partition), val(chrom), val(cellType), path(hic_out), path(hic_cell), path(unmatched_rows)

//     output:
//     tuple val(identifier), val(partition), val(cellType), path("*.csv"), val(chrom), path("${hic_input_name}.${cellType}.${partition}.README.txt")

//     script:
//     """
//     hic_out_list=($hic_out)
//     IFS=\$'\n' hic_out_list=(\$(sort -V <<<"\${hic_out_list[*]}")); unset IFS
//     outf=`echo "\${hic_out_list[0]}" | awk -F. '\$1=\$1{for(i=1; i<=(NF-4); ++i) printf \$i""FS; print "" }'`
//     ext=`echo "\${hic_out_list[0]}" | awk -F. '\$1=\$1{for(i=(NF-2); i<=(NF-1); ++i) printf \$i""FS; print "" }'`
//     outf="\${outf}\${ext}${partition}.csv"
//     for file in "\${hic_out_list[@]}"
//     do
//         tail -n +2 \$file >> tmp_all.csv
//     done
//     head -n1 \${hic_out_list[0]} > \$outf
//     cat tmp_all.csv >> \$outf

//     hic_out_list=($hic_cell)
//     IFS=\$'\n' hic_out_list=(\$(sort -V <<<"\${hic_out_list[*]}")); unset IFS
//     outf=`echo "\${hic_out_list[0]}" | awk -F. '\$1=\$1{for(i=1; i<=(NF-2); ++i) printf \$i""FS; print "" }'`
//     outf="\${outf}pos.${partition}.csv"
//     for file in "\${hic_out_list[@]}"
//     do
//         tail -n +2 \$file >> tmp.csv
//     done
//     head -n1 \${hic_out_list[0]} > \$outf
//     cat tmp.csv >> \$outf

//     hic_out_list=($unmatched_rows)
//     IFS=\$'\n' hic_out_list=(\$(sort -V <<<"\${hic_out_list[*]}")); unset IFS
//     outf=`echo "\${hic_out_list[0]}" | awk -F. '\$1=\$1{for(i=1; i<=(NF-2); ++i) printf \$i""FS; print "" }'`
//     outf="\${outf}pos.${partition}.csv"
//     for file in "\${hic_out_list[@]}"
//     do
//         tail -n +2 \$file >> tmp_unmat.csv
//     done
//     head -n1 \${hic_out_list[0]} > \$outf
//     cat tmp_unmat.csv >> \$outf

//     echo "$identifier" > ${hic_input_name}.${cellType}.${partition}.README.txt

//     rm tmp_all.csv tmp.csv tmp_unmat.csv
//     """
// }

process COMBINE_FILES{
    publishDir "${params.publish_dir}"

    input:
    tuple val(identifier), path(hic_out), val(outf)

    output:
    tuple val(identifier), path("${outf}")

    script:
    """
    hic_out_list=($hic_out)
    hic_out_list=(\$(sort -V <<<"\${hic_out_list[*]}"))
    for file in "\${hic_out_list[@]}"
    do
        tail -n +2 \$file >> tmp.csv
    done
    head -n1 \${hic_out_list[0]} > $outf
    cat tmp.csv >> $outf
    """
}

process RANDOMIZE {
    stageInMode 'link'
    label 'process_intense'

    input:
    tuple val(cellType), val(label), path(data_h5), val(chunk), val(n_chunks)

    output:
    tuple val(cellType), val(label), path("data.${cellType}.rep${label}.rnd.${chunk}.h5")

    script:
    """
    #!/usr/bin/env python

    import math
    import h5py
    import pandas as pd
    import numpy as np

    with h5py.File("${data_h5}",'r') as h5f0:
        dset_size = h5f0['label'].shape[0]
        NUM_SEQ = h5f0['enh_seq'].shape[2]
        NUM_REP = h5f0['enh_dnase'].shape[2]

        np.random.seed(42)
        permut = np.random.permutation(dset_size)
        permut_inv = np.array(sorted(range(len(permut)), key=lambda k: permut[k]))

        n_chunks = $n_chunks
        chunk = $chunk
        assert chunk<n_chunks

        init_index = chunk*math.floor(dset_size/n_chunks) + min(chunk,dset_size%n_chunks)
        last_index = init_index + math.floor(dset_size/n_chunks) + int(chunk<(dset_size%n_chunks))

        indexes = permut_inv[init_index:last_index]
        indexes_sort = np.array(sorted(range(len(indexes)), key=lambda k: indexes[k]))
        indexes_sort_inv = np.array(sorted(range(len(indexes)), key=lambda k: indexes_sort[k]))

        with h5py.File(f"data.${cellType}.rep${label}.rnd.{chunk}.h5",'w') as h5f:
            for k in h5f0.keys():
                h5f.create_dataset(k, data=h5f0[k][indexes[indexes_sort]][indexes_sort_inv])
            # for k in h5f0.keys():
            #     h5f.create_dataset(k, (len(indexes),)+h5f0[k].shape[1:], dtype = h5f0[k].dtype)
            # for ind, index in enumerate(indexes):
            #     if ind%100==0:
            #         print(f"{ind}/{len(indexes)}", flush=True)
            #     for k in h5f0.keys():
            #         h5f[k][ind] = h5f0[k][index]
    """
}


process KNAPSACK{
    input:
    tuple val(identifier), val(list_count), val(frac)

    output:
    tuple val(identifier), path("select_indices.csv")

    script:    
    """
    knapsack.py '$list_count' '$frac'
    # knapsack.py '[37101, 79501, 122501, 77100, 119001, 66701]' '[0.2,0.2,0.6]'
    """
}


process CHUNK_DATA {

    input:
    tuple val(cellType), val(label), path(hic_aug), path(enh_index), path(pr_index), val(n_chunks)

    output:
    tuple val(cellType), val(label), path("${hic_aug}.*"), path("${enh_index}.*"), path("${pr_index}.*")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    df=pd.read_csv("${hic_aug}")
    df_enhindex = pd.read_csv("${enh_index}",header=None)
    df_prindex = pd.read_csv("${pr_index}",header=None)

    N=df.shape[0]
    chunks=${n_chunks}
    N_=np.array([N//chunks]*chunks)
    N_[:N%chunks] += 1


    prev_rows=0
    for i, r in enumerate(N_):
        df.iloc[prev_rows:prev_rows+r].to_csv(f"${hic_aug}.{i}", index=False)
        df_enhindex.iloc[prev_rows:prev_rows+r].to_csv(f"${enh_index}.{i}", index=False,header=None)
        df_prindex.iloc[prev_rows:prev_rows+r].to_csv(f"${pr_index}.{i}", index=False,header=None)
        prev_rows += r
    """
}

process SUBLIST{
    input:
    tuple val(identifier), val(select_list), val(list_index)

    output:
    tuple val(identifier) , val(selected_list)

    script:
    _list_index = list_index.tokenize(', []').collect{cc -> cc.toInteger()}
    // println "${select_list}"
    _select_list = select_list.toList()
    selected_list = _list_index.collect{ x -> _select_list[x]}
    left_elems = _select_list.minus(selected_list)
    """    
    # echo '$identifier'
    # echo '$select_list'
    echo '$_list_index'

    # echo 'out:'
    # echo "selected_list"
    # echo "left_elems"
    """
}


// DNase/ATAC preprocessing step 3.1
// generate chromatin openness core profile (COscore) for promoter
process CHROM_OPENN_SCORE_PROFILE {
    publishDir "${params.publish_dir}"

    input:
    tuple val(cell), val(rep), path(npzlist), path(regelem), path(regelem_bg), val(regelem_name)

    output:
    tuple val(cell), val(rep), path("${regelem_name}_COscore.${cell}.rep${rep}.h5.gz")

    script:  
    // npzlist_ = npzlist.toList()
    """
    python -c "import preptools as pt; \
    pt.generate_dnase('$npzlist'.split(), \
        '$regelem', \
        '$regelem_bg', \
        $params.bgWindow, \
        '${regelem_name}_COscore.${cell}.rep${rep}')"
    gzip ${regelem_name}_COscore.${cell}.rep${rep}.h5
    """
}
