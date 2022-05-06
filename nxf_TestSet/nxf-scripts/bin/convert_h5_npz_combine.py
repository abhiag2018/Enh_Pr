#!/usr/bin/env python

import sys

import glob
import natsort
import math
import h5py
import pandas as pd
import numpy as np

n_chunks = 100

cell_list = sys.argv[1:]
cell='-'.join(cell_list)
print(cell, flush=True)
cell_data=[]
for c in cell_list:
    cell_data += glob.glob(f"{c}/data_*/d*_*.npz")

cell_data=np.array(cell_data)

RESIZED_LEN=2000
NUM_REP=3
NUM_SEQ=4

with h5py.File(f"{cell}/{cell}.h5",'w') as h5f0:
    enh_dnase = h5f0.create_dataset('enh_dnase', (len(cell_data),1,NUM_REP,RESIZED_LEN), dtype='float64')    
    pr_dnase = h5f0.create_dataset('pr_dnase', (len(cell_data),1,NUM_REP,1000), dtype='float64')    
    enh_seq = h5f0.create_dataset('enh_seq', (len(cell_data),1,NUM_SEQ,RESIZED_LEN), dtype='float64')    
    pr_seq = h5f0.create_dataset('pr_seq', (len(cell_data),1,NUM_SEQ,1000), dtype='float64')    
    label = h5f0.create_dataset('label', (len(cell_data),), dtype='float64')    

    for chunk in range(n_chunks):
        init_index = chunk*math.floor(len(cell_data)/n_chunks) + min(chunk,len(cell_data)%n_chunks)
        last_index = init_index + math.floor(len(cell_data)/n_chunks) + int(chunk<(len(cell_data)%n_chunks))

        with h5py.File(f"{cell}/{cell}.{chunk}.h5",'r') as h5f:
            # h5f0['enh_dnase'][init_index:last_index] = data_np['enh_dnase']
            enh_dnase.write_direct(h5f['enh_dnase'][:], dest_sel=np.s_[init_index:last_index])
            pr_dnase.write_direct(h5f['pr_dnase'][:], dest_sel=np.s_[init_index:last_index])
            enh_seq.write_direct(h5f['enh_seq'][:], dest_sel=np.s_[init_index:last_index])
            pr_seq.write_direct(h5f['pr_seq'][:], dest_sel=np.s_[init_index:last_index])
            label.write_direct(h5f['label'][:], dest_sel=np.s_[init_index:last_index])
