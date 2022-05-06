#!/usr/bin/env python

import h5py
import pandas as pd
import numpy as np
import sys
import glob

dpath=f"/projects/li-lab/agarwa/CUBE/DeepTact/DeepTact-Pipeline/storeDir/{sys.argv[1]}"
RESIZED_LEN=2000
NUM_REP=3
NUM_SEQ=4
D=np.sum([pd.read_csv(f).shape[0] for f in glob.glob(f"{dpath}.*.csv")])
indexes = np.random.permutation(D)

with h5py.File(f"{dpath}.h5",'w') as h5f:
    h5f.create_dataset('enh_dnase', (D,1,NUM_REP,RESIZED_LEN), dtype='float64')    
    h5f.create_dataset('pr_dnase', (D,1,NUM_REP,1000), dtype='float64')    
    h5f.create_dataset('enh_seq', (D,1,NUM_SEQ,RESIZED_LEN), dtype='float64')    
    h5f.create_dataset('pr_seq', (D,1,NUM_SEQ,1000), dtype='float64')    
    h5f.create_dataset('label', (D,), dtype='float64')    

    ind=0
    for chunk in range(len(glob.glob(f"{dpath}.*.csv"))):
        DF = pd.read_csv(f"{dpath}.{chunk}.csv")
        DF[DF.columns.drop('label')] = DF[DF.columns.drop('label')].applymap(eval).applymap(list)

        DF['enh_dnase'] = DF.apply(lambda row: np.array([row[f'enh_dnase_{x}'] for x in range(NUM_REP)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)
        DF['pr_dnase'] = DF.apply(lambda row: np.array([row[f'pr_dnase_{x}'] for x in range(NUM_REP)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)
        DF['enh_seq'] = DF.apply(lambda row: np.array([row[f'enh_seq_{x}'] for x in range(NUM_SEQ)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)
        DF['pr_seq'] = DF.apply(lambda row: np.array([row[f'pr_seq_{x}'] for x in range(NUM_SEQ)], dtype=np.dtype('float64'))[np.newaxis,:,:],axis=1)

        DF = DF[['enh_dnase', 'pr_dnase', 'enh_seq', 'pr_seq','label']]

        for _, row in DF.iterrows():
            h5f['enh_dnase'][indexes[ind]] = row['enh_dnase']
            h5f['pr_dnase'][indexes[ind]] = row['pr_dnase']
            h5f['enh_seq'][indexes[ind]] = row['enh_seq']
            h5f['pr_seq'][indexes[ind]] = row['pr_seq']
            h5f['label'][indexes[ind]] = row['label']
            ind += 1

pd.DataFrame({'indexes':indexes}).to_csv(f"{dpath}-order.csv", index=False)
# ## Read

# with h5py.File(f"{ddir}/dataset-train.h5",'r') as h5f:
#     for i in h5f['Label'].shape[0]:
#         yield {'enh_dnase': h5f['enh_dnase'][i], 'pr_dnase': h5f['pr_dnase'][i], 'enh_seq': h5f['enh_seq'][i], 'pr_seq': h5f['pr_seq'][i]}, h5f['Label'][i]