#!/usr/bin/env python

import numpy as np
import sklearn
import sklearn.metrics
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_roc(saveto, name, labels, predictions):
    mpl.rcParams['figure.figsize'] = (12, 10)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fp, tp, _ = sklearn.metrics.roc_curve(labels, predictions)

    auc = sklearn.metrics.auc(fp, tp)
    print("auroc:", auc)
    plt.plot(100*fp, 100*tp, label = name+f"; auc : {auc:.2f}", linewidth=2, color=colors[0], linestyle='--')
    plt.legend()
    plt.xlabel('False positives [%]')
    plt.ylabel('True positives [%]')
    plt.xlim([0,100])
    plt.ylim([0,100])
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(saveto)
    plt.close('all')

def plot_prc(saveto, name, labels, predictions):
    mpl.rcParams['figure.figsize'] = (12, 10)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    pr, rec, _ = sklearn.metrics.precision_recall_curve(labels, predictions)

    auc = sklearn.metrics.auc(rec, pr)
    print("auprc:", auc)
    plt.plot(100*rec, 100*pr, label = name+f"; auc : {auc:.2f}", linewidth=2, color=colors[0], linestyle='--')
    plt.legend()
    plt.xlabel('Recall [%]')
    plt.ylabel('Precision [%]')
    plt.xlim([0,100])
    plt.ylim([0,100])
    plt.grid(True)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.savefig(saveto)
    plt.close('all')


def eval_stats(eval_npz, cell="type-1", eval_cell="type-2",out_prefix=""):
    # cell=CELL.split('/')[-1]
    legend_str = f"train = {cell}; test = {eval_cell}"

    pred = np.load(eval_npz)['pred_score']
    gt = np.load(eval_npz)['true_labels']

    plot_roc(f"{out_prefix}_roc.png", legend_str, gt, pred)
    plot_prc(f"{out_prefix}_prc.png", legend_str, gt, pred)

    print("f1-score", sklearn.metrics.f1_score(gt, pred>=0.5))
    print("acc", sklearn.metrics.accuracy_score(gt, pred>=0.5))

    C = sklearn.metrics.confusion_matrix(gt, pred>=0.5)
    num_dig = np.ceil(np.log10(max(C.flatten())))
    print()
    print(f"TN={C[0,0]<num_dig} FP={C[0,1]<num_dig} ")
    print(f"FN={C[1,0]<num_dig} TP={C[1,1]<num_dig} ")
