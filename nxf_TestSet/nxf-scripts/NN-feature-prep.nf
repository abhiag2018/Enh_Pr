nextflow.preview.dsl=2
include { CONSOLIDATE_OUT as CONSOLIDATE_1 } from './module'
include { CONSOLIDATE_OUT as CONSOLIDATE_2 } from './module'
include { RANDOMIZE as RANDOMIZE_1 } from './module'
include { RANDOMIZE as RANDOMIZE_2 } from './module'
include { RANDOMIZE as RANDOMIZE_3 } from './module'
include { RANDOMIZE as RANDOMIZE_4 } from './module'
include { RANDOMIZE as RANDOMIZE_A } from './module'

include { CHUNK_DATA as CHUNK_DATA_1 } from './module'
include { CHUNK_DATA as CHUNK_DATA_2 } from './module'
include { CHUNK_DATA as CHUNK_DATA_3 } from './module'
include { CHUNK_DATA as CHUNK_DATA_4 } from './module'
include { CHUNK_DATA as CHUNK_DATA_A } from './module'


Channel.fromPath("$params.dnaseq_inp")
  .splitCsv(header:true, sep:',')
  .map { row -> [ file("$params.store_dir/$row.enhancer", checkIfExists: true), file("$params.store_dir/$row.promoter", checkIfExists: true) ]  }
  .set {ch_dnaseq}

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

Channel.fromPath("$params.coscore_inp")
  .splitCsv(header:true, sep:',')
  .map { row -> [ row.celltype, row.repetition, row.regelem, file("$row.h5", checkIfExists: true) ]  }
  .set {ch_co_score}

ch_enh_co_score=ch_co_score.filter{ it[2] == "enhancer" }.map{ it -> [it[0],it[1],it[3]]}
ch_pr_co_score=ch_co_score.filter{ it[2] == "promoter" }.map{ it -> [it[0],it[1],it[3]]}

Channel.fromPath("$params.pchic_inp")
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.celltype, row.label, file("$row.hic", checkIfExists: true) ]  }
    .set { ch_hic_aug }
if (params.dev){
    ch_hic_aug = ch_hic_aug
        .splitText(by:1000,keepHeader:true,file:true)
        .take(1)
}

// align enhancers and promoters in the pchic interaction data with the promoter and enhancer lists
process ALIGN_PCHIC_ENHPR {

    input:
    tuple val(cellType), val(label), path(hic_aug), path(enhancer_DNAseq), path(promoter_DNAseq)

    output:
    tuple val(cellType), val(label), path(hic_aug), path("enh.index"), path("pr.index") 

    script:  
    """
    python -c "import pandas as pd, numpy as np; pd.DataFrame({'name':np.load('$enhancer_DNAseq',allow_pickle=True)['name']}).to_csv('enh_name.csv')"
    python -c "import pandas as pd, numpy as np; pd.DataFrame({'name':np.load('$promoter_DNAseq',allow_pickle=True)['name']}).to_csv('pr_name.csv')"
    awk -F, 'NR==FNR{a[\$2]=\$1; next} {if (\$4 in a) print a[\$4]}' enh_name.csv $hic_aug > enh.index
    awk -F, 'NR==FNR{a[\$2]=\$1; next} {if (\$8 in a) print a[\$8]}' pr_name.csv $hic_aug > pr.index
    """
}


process COMBINE_PCHIC_FEATURES {
    tag "data.${cellType}.{i}.tar.gz"
    label 'process_intense'

    input:
    tuple val(cellType), val(label), path(hic_aug), path(enh_index), path(pr_index), path(enh_DNAseq), path(pr_DNAseq), path(enh_coscore_gz), path(pr_coscore_gz)

    output:
    // tuple val(cellType), path("data.${cellType}.*.h5")
    tuple val(cellType), val(label), path("data.${cellType}.*.h5")

    script:
    index=pr_index.getExtension()
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    import subprocess
    import h5py
    import os

    pos_data = np.array([],dtype=str)
    neg_data = np.array([],dtype=str)
    enh_coscore = "$enh_coscore_gz".split()
    pr_coscore = "$pr_coscore_gz".split()

    if (enh_coscore[0].split('.')[-1]=="gz"):
        if os.system("gzip -df $enh_coscore_gz $pr_coscore_gz") != 0:
            raise Exception("cannot de-compress")
        else:
            enh_coscore = ['.'.join(x.split('.')[:-1]) for x in enh_coscore]
            pr_coscore = ['.'.join(x.split('.')[:-1]) for x in pr_coscore]

    hicInteractions = pd.read_csv('$hic_aug')
    enh_index =  np.array([int(x) for x in open("$enh_index", "r").readlines()])
    pr_index = np.array([int(x) for x in open("$pr_index", "r").readlines()])
    enh_sort = np.array(sorted(range(len(enh_index)), key=lambda k: enh_index[k]))
    enh_sortinv = np.array(sorted(range(len(enh_index)), key=lambda k: enh_sort[k]))
    pr_sort = np.array(sorted(range(len(pr_index)), key=lambda k: pr_index[k]))
    pr_sortinv = np.array(sorted(range(len(pr_index)), key=lambda k: pr_sort[k]))

    enh_dnaseq=np.load("$enh_DNAseq",allow_pickle=True)
    pr_dnaseq=np.load("$pr_DNAseq",allow_pickle=True)
    enh_reps = []
    for enh_rep in enh_coscore:
        enh_reps.append(np.array(h5py.File(enh_rep, 'r')['expr'])[enh_index])
    pr_reps = []
    for pr_rep in pr_coscore:
        pr_reps.append(np.array(h5py.File(pr_rep, 'r')['expr'])[pr_index])

    assert len(pr_index) == len(enh_index)
    assert len(hicInteractions) == len(enh_index)
    assert np.sum(hicInteractions.label.apply(lambda x:not(x==1 or x==0))) == 0

    assert enh_dnaseq['sequence'].shape[0] == h5py.File(enh_coscore[0],'r')['expr'].shape[0]
    assert pr_dnaseq['sequence'].shape[0] == h5py.File(pr_coscore[0],'r')['expr'].shape[0]

    _enh_dnaseq = enh_dnaseq['sequence']
    _pr_dnaseq = pr_dnaseq['sequence']
    _enh_reps = np.stack(enh_reps,axis=1)
    _pr_reps = np.stack(pr_reps,axis=1)
    _pr_loc = pr_dnaseq['loc']
    _pr_name = pr_dnaseq['name']
    _enh_loc = enh_dnaseq['loc']
    _enh_name = enh_dnaseq['name']

    NUM_REP = _enh_reps.shape[1]
    NUM_SEQ = 4

    enh_start = np.array(list(map(lambda x:int(x[0]),map(lambda x:x.split('-'),map(lambda x:x.split(':')[1], _enh_loc[enh_index])))))
    pr_start = np.array(list(map(lambda x:int(x[0]),map(lambda x:x.split('-'),map(lambda x:x.split(':')[1], _pr_loc[pr_index])))))

    enh_shift = np.array(hicInteractions.enhancer_start - enh_start)
    pr_shift = np.array(hicInteractions.promoter_start - pr_start)

    def slice(X,win_shift,len_win):
        return X[np.arange(X.shape[0])[:,None,None,None], np.arange(X.shape[1])[None,:,None,None], np.arange(X.shape[2])[None,None,:,None], np.add.outer(enh_shift,np.arange(len_win))[:,None,None,:] ]

    # create hdf5 chunks
    with h5py.File(f'data.${cellType}.${index}.h5', 'w') as h5f:
        h5f.create_dataset('enh_seq', data=slice(_enh_dnaseq[enh_index[enh_sort]][enh_sortinv].reshape((len(enh_index),1,-1,4)).transpose(0,1,3,2),enh_shift,$params.enhancer_window) , dtype='float64')
        h5f.create_dataset('pr_seq', data=slice(_pr_dnaseq[pr_index[pr_sort]][pr_sortinv].reshape((len(enh_index),1,-1,4)).transpose(0,1,3,2),pr_shift,$params.promoter_window) , dtype='float64')
        h5f.create_dataset('enh_dnase', data=slice(_enh_reps[:,np.newaxis,:,:],enh_shift, $params.enhancer_window)  , dtype='float64')
        h5f.create_dataset('pr_dnase', data=slice(_pr_reps[:,np.newaxis,:,:],pr_shift,$params.promoter_window)  , dtype='float64')
        h5f.create_dataset('label', data=hicInteractions['label'], dtype='float64')
        h5f.create_dataset('enh.loc', data=_enh_loc[enh_index[enh_sort]][enh_sortinv], dtype=h5py.special_dtype(vlen=str))
        h5f.create_dataset('pr.loc', data=_pr_loc[pr_index[pr_sort]][pr_sortinv], dtype=h5py.special_dtype(vlen=str))
        h5f.create_dataset('enh.name', data=_enh_name[enh_index[enh_sort]][enh_sortinv], dtype=h5py.special_dtype(vlen=str))
        h5f.create_dataset('pr.name', data=_pr_name[pr_index[pr_sort]][pr_sortinv], dtype=h5py.special_dtype(vlen=str))

        # h5f.create_dataset('enh_dnase', (len(enh_index),1,NUM_REP,$params.enhancer_window), dtype='float64')    
        # h5f.create_dataset('pr_dnase', (len(enh_index),1,NUM_REP,$params.promoter_window), dtype='float64')    
        # h5f.create_dataset('enh_seq', (len(enh_index),1,NUM_SEQ,$params.enhancer_window), dtype='float64')    
        # h5f.create_dataset('pr_seq', (len(enh_index),1,NUM_SEQ,$params.promoter_window), dtype='float64')    
        # h5f.create_dataset('label', (len(enh_index),), dtype='float64')    
        # h5f.create_dataset('enh.loc', (len(enh_index),),  dtype=h5py.special_dtype(vlen=str))
        # h5f.create_dataset('enh.name', (len(enh_index),),  dtype=h5py.special_dtype(vlen=str))
        # h5f.create_dataset('pr.loc', (len(enh_index),),  dtype=h5py.special_dtype(vlen=str))
        # h5f.create_dataset('pr.name', (len(enh_index),),  dtype=h5py.special_dtype(vlen=str))

        # for i in range(len(enh_index)):
        #     if i %100 ==0:
        #         print(i, flush=True)
        #     enh_start, _ = map(int, _enh_loc[enh_index[i]].split(':')[1].split('-'))
        #     pr_start, _ = map(int, _pr_loc[pr_index[i]].split(':')[1].split('-'))

        #     enh_shift = hicInteractions.enhancer_start[i] - enh_start
        #     pr_shift = hicInteractions.promoter_start[i] - pr_start

        #     h5f['enh_seq'][i] = _enh_dnaseq[enh_index[i]].reshape((1,-1,4)).transpose(0,2,1)[:,:,enh_shift:enh_shift+$params.enhancer_window] 
        #     h5f['pr_seq'][i] = _pr_dnaseq[pr_index[i]].reshape((1,-1,4)).transpose(0,2,1)[:,:,pr_shift:pr_shift+$params.promoter_window] 
        #     h5f['enh_dnase'][i] = _enh_reps[i][np.newaxis,:,:][:,:,enh_shift:enh_shift+$params.enhancer_window] 
        #     h5f['pr_dnase'][i] = _pr_reps[i][np.newaxis,:,:][:,:,pr_shift:pr_shift+$params.promoter_window] 
        #     h5f['label'][i] = hicInteractions['label'][i]
        #     h5f['enh.loc'][i] = _enh_loc[enh_index[i]]
        #     h5f['pr.loc'][i] = _pr_loc[pr_index[i]]
        #     h5f['enh.name'][i] = _enh_name[enh_index[i]]
        #     h5f['pr.name'][i] = _pr_name[pr_index[i]]
    """
}

workflow combine_data_v1{
    take:
    ch_enh_co_score
    ch_enhancer_bed_prep
    ch_pr_co_score
    ch_promoter_bed_prep
    ch_dnaseq
    ch_hic_aug_0

    main:
    ch_coscore = ch_enh_co_score.join(ch_pr_co_score, by:[0,1]).groupTuple(by:0).map(it -> [it[0], it[2], it[3]])

    ch_hic_aug = ch_hic_aug_0.filter{ file(it[2]).countLines()>1 }

    ch_index = ALIGN_PCHIC_ENHPR(ch_hic_aug.combine(ch_dnaseq))
    ch_index_ = CHUNK_DATA_1(ch_index.filter{ it[1] =~ ~/train$/ }.combine(Channel.value(params.chunks_train)))
        .mix(CHUNK_DATA_2(ch_index.filter{ it[1] =~ ~/test$/ }.combine(Channel.value(params.chunks_test))))
        .mix(CHUNK_DATA_3(ch_index.filter{ it[1] =~ ~/val$/ }.combine(Channel.value(params.chunks_val))))
        .mix(CHUNK_DATA_4(ch_index.filter{ it[1] =~ ~/left-overs$/ }.combine(Channel.value(params.chunks_val))))
        .mix(CHUNK_DATA_A(ch_index.filter{ it[1] =~ ~/all$/ }.combine(Channel.value(params.chunks_all))))


    ch_comb = ch_index_.transpose().combine(ch_dnaseq).combine(ch_coscore,by:0)

    // ch_comb = ch_comb.filter{ it[0] =~ "type93" }.filter{ it[3].getExtension() =~ "185|197|216" }
    // ch_comb = ch_comb.filter{ it[0] =~ "type85" }.filter{ it[3].getExtension() =~ "208" }
    // ch_comb = ch_comb.filter{ it[0] =~ "type83" }.filter{ it[3].getExtension() =~ "167" }

    ch_deept_features = COMBINE_PCHIC_FEATURES(ch_comb).groupTuple(by:[0,1]).combine(Channel.value("")) | CONSOLIDATE_1

    if (params.randomize){
        chunks_leftovers = params.chunks_val
        ch_deept_features_random = RANDOMIZE_1(ch_deept_features.filter{ it[1] =~ ~/train$/ }.combine(Channel.fromList(0..(params.chunks_train-1))).combine(Channel.value(params.chunks_train)))
            .mix(RANDOMIZE_2(ch_deept_features.filter{ it[1] =~ ~/test$/ }.combine(Channel.fromList(0..(params.chunks_test-1))).combine(Channel.value(params.chunks_test))))
            .mix(RANDOMIZE_3(ch_deept_features.filter{ it[1] =~ ~/val$/ }.combine(Channel.fromList( 0..(params.chunks_val-1))).combine(Channel.value(params.chunks_val))))
            .mix(RANDOMIZE_4(ch_deept_features.filter{ it[1] =~ ~/left-overs$/ }.combine(Channel.fromList( 0..(chunks_leftovers-1))).combine(Channel.value(chunks_leftovers))))
            .mix(RANDOMIZE_A(ch_deept_features.filter{ it[1] =~ ~/all$/ }.combine(Channel.fromList( 0..(params.chunks_all-1))).combine(Channel.value(params.chunks_all))))
            .groupTuple(by:[0,1])
            
        ch_deept_features_random = CONSOLIDATE_2(ch_deept_features_random.combine(Channel.value(".rnd")))
    }else{
        ch_deept_features_random = ch_deept_features
    }

    emit: ch_deept_features_random
}

workflow{
    combine_data_v1(ch_enh_co_score, ch_enhancer_bed_prep, ch_pr_co_score, ch_promoter_bed_prep, ch_dnaseq, ch_hic_aug).view()
}













