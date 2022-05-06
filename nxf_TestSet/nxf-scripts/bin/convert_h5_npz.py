#!/usr/bin/env python

import sys

import glob
import natsort
import math
import h5py
import pandas as pd
import numpy as np
np.random.seed(42)

chunk = int(sys.argv[1])
n_chunks = 100
assert chunk<n_chunks

cell_list = sys.argv[2:]
cell='-'.join(cell_list)
print(cell, flush=True)
cell_data=[]
for c in cell_list:
    cell_data += glob.glob(f"{c}/data_*/d*_*.npz")

cell_data=np.array(natsort.natsorted(cell_data))

RESIZED_LEN=2000
NUM_REP=3
NUM_SEQ=4

permut = np.random.permutation(len(cell_data))
permut_inv = np.array(sorted(range(len(permut)), key=lambda k: permut[k]))


init_index = chunk*math.floor(len(cell_data)/n_chunks) + min(chunk,len(cell_data)%n_chunks)
last_index = init_index + math.floor(len(cell_data)/n_chunks) + int(chunk<(len(cell_data)%n_chunks))

indexes = permut_inv[init_index:last_index]
# permut[indexes]



with h5py.File(f"{cell}.{chunk}.h5",'w') as h5f:
    h5f.create_dataset('enh_dnase', (len(indexes),1,NUM_REP,RESIZED_LEN), dtype='float64')    
    h5f.create_dataset('pr_dnase', (len(indexes),1,NUM_REP,1000), dtype='float64')    
    h5f.create_dataset('enh_seq', (len(indexes),1,NUM_SEQ,RESIZED_LEN), dtype='float64')    
    h5f.create_dataset('pr_seq', (len(indexes),1,NUM_SEQ,1000), dtype='float64')    
    h5f.create_dataset('label', (len(indexes),), dtype='float64')    

    for ind, data in enumerate(cell_data[indexes]):
        if ind%100==0:
            print(f"{ind}/{len(indexes)}", flush=True)
        data_np = np.load(data)
        h5f['enh_dnase'][ind] = data_np['enh_dnase']
        h5f['pr_dnase'][ind] = data_np['pr_dnase']
        h5f['enh_seq'][ind] = data_np['enh_seq']
        h5f['pr_seq'][ind] = data_np['pr_seq']
        h5f['label'][ind] = data_np['label']

if chunk==0:
    pd.DataFrame({'indexes':cell_data[permut]}).to_csv(f"{cell}-order.csv", index=False)
# ## Read

# with h5py.File(f"{ddir}/dataset-train.h5",'r') as h5f:
#     for i in h5f['Label'].shape[0]:
#         yield {'enh_dnase': h5f['enh_dnase'][i], 'pr_dnase': h5f['pr_dnase'][i], 'enh_seq': h5f['enh_seq'][i], 'pr_seq': h5f['pr_seq'][i]}, h5f['Label'][i]