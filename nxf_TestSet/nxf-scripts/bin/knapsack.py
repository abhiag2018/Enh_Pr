#!/usr/bin/env python
import sys
from multiset import *
import numpy as np
import pandas as pd

list_count = eval(sys.argv[1])
frac = eval(sys.argv[2])

def knapSack(W, wt): 
    print("W, wt:", W, np.array(wt), flush=True)

    val = wt
    n = len(wt)

    K = np.zeros((n+1,W+1))#[[0 for x in range(W + 1)] for x in range(n + 1)] 
    Index = np.zeros((n+1,W+1))
    Wt_Index = np.zeros((n+1,W+1))

    # Build table K[][] in bottom up manner 
    for i in range(n + 1): 
        for w in range(W + 1): 
            if i == 0 or w == 0: 
                Index[i][w] = None
                Wt_Index[i][w] = -1
                K[i][w] = 0
            elif wt[i-1] <= w: 
                if val[i-1] + K[i-1][w-wt[i-1]] >=  K[i-1][w]:
                    Index[i][w] = i-1
                    Wt_Index[i][w] = w-wt[i-1]
                    K[i][w] = val[i-1] + K[i-1][w-wt[i-1]]
                else:
                    Index[i][w] = None
                    Wt_Index[i][w] = w
                    K[i][w] = K[i-1][w]
            else: 
                Index[i][w] = None
                Wt_Index[i][w] = w
                K[i][w] = K[i-1][w] 

    list_items = []
    wt_index = 0
    ind = n; wt_index = W
    while (wt_index >= 0 and ind>=0):
        # print(ind,wt_index)
        _ind  = Index[ind,int(wt_index)]
        _wt_index = Wt_Index[ind,int(wt_index)]
        list_items = list_items + ([] if np.isnan(_ind) else [_ind])
        ind = ind-1
        wt_index = _wt_index

    list_items = [int(x) for x in list_items]
    wt = np.array(wt)
    # print('items')
    # print(wt)
    # print(wt[list_items])
    # print(sum(wt[list_items]),  K[n][W])

    return list_items, K[n][W]

wt = np.array([x-1 for x in list_count])
wt_sum = np.sum(wt)
# wt = np.array([x-1 for x in [37101, 79501, 122501, 77100, 119001, 66701]])#, 80400, 94601, 109401, 56901, 119601, 1, 49001, 205601, 13201, 99301, 183201, 124501, 60501, 125801, 45901]])
list_ind_orig = np.arange(len(wt))

fractions = {'test':frac[0], 'val':frac[1] ,'train':frac[2]}
print(fractions, flush=True)

assert np.sum(list(fractions.values()))==1
list_rep_ind_ = {}
W_ = {}
while len(fractions)>0:
    _key = max(fractions, key=fractions.get)

    W = int(fractions[_key]*wt_sum)
    list_items, W = knapSack(W, wt)

    W_[_key] = int(W)
    list_rep_ind_[_key] = list_ind_orig[list_items]
    list_ind_orig = np.array([list_ind_orig[k] for k in (set(range(len(wt))) - set(list_items))])
    wt = [wt[k] for k in (set(range(len(wt))) - set(list_items))]

    print(W_, list_rep_ind_, flush=True,end="\n\n\n")
    fractions.pop(_key)

print("left-overs:", list_ind_orig)
W_['left-overs'] = -1
list_rep_ind_['left-overs'] = list_ind_orig

labels = ['test','val','train','left-overs']
df = pd.DataFrame({'list_index':map(lambda x:list(list_rep_ind_[x]), labels),'num_cols':map(lambda x:W_[x], labels),'label':labels})
df.to_csv("select_indices.csv",sep=';',index=False)


