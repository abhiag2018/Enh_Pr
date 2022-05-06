#!/usr/bin/env python
import colored_traceback.always
import traceback
import re
import bisect

import numpy as np
import h5py
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

import gtfparse
import deeptools.countReadsPerBin as crpb

def getBamCounts(bam_file, chrom, chromLen, outputf = None):
    """
    read a bam file and return the readCounts for position - chrom:0-chromLen
    inputs: 
        chrom : chromosome to read
        chromLen : scalar int
    """
    cr = crpb.CountReadsPerBin([bam_file], binLength=1, stepSize=1)
    try:
        arr = cr.count_reads_in_region(chrom[3:], 0, chromLen)[0]
    except NameError:
        arr = cr.count_reads_in_region(chrom, 0, chromLen)[0]

    if outputf:
        np.savez(outputf, count=arr)
        # np.load(outputf,allow_pickle=True)['count']
    return arr


def generate_dnase(chr_list, bed_path, bg_path, bgWin, output_fn):
    # chr_list = chrom_npz_fileList.split()

    chrData = {}
    for chr_npz in chr_list:
        chrom = chr_npz.split(".")[-2]
        chrData[chrom] = np.load(chr_npz,allow_pickle=True)['count']

    prDF = pd.read_csv(bed_path, delimiter='\t', header=None)
    prDF_bg = pd.read_csv(bg_path, delimiter='\t', header=None)
    
    dnase_profile = np.array(prDF.apply(lambda df:chrData[df[0]][df[1]:df[2]].flatten(), axis=1).tolist())
    bg_profile = np.array(prDF_bg.apply(lambda df:np.sum(chrData[df[0]][df[1]:df[2]]),axis=1).tolist())

    _dnase_profile = np.sum(dnase_profile,axis=1)
    assert np.sum(_dnase_profile[bg_profile==0])==0
    assert np.logical_not(np.logical_xor((bg_profile<1), (bg_profile==0))).all()
    bg_profile[bg_profile==0]=1
    out = dnase_profile / bg_profile[:,None] * bgWin

    with h5py.File(f"{output_fn}.h5", 'w') as hf:
        hf.create_dataset('expr',  data=out)
    # out = h5py.File(f"{output_fn}.h5",'r')['expr']

    # np.savez(f"{output_fn}.npz",expr = out)
    # np.load(f"{output_fn}.npz",allow_pickle=True)['expr']
    return 0


def process_promoter_v1(gtf_input, window, all_chroms, out_fname):
    """
    process gtf files for hg19 gtf ensembl release 75 fiels
    """
    gtf = gtfparse.read_gtf(gtf_input)
    gtf['seqname'] = gtf['seqname'].apply(lambda x:"chr"+x)
    gtf = gtf[gtf['seqname'].apply(lambda x:x in all_chroms)]
    cols_select = ['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'transcript_id', 'transcript_name', 'protein_id', 'score']
    gtf1 = gtf[(gtf.feature == "transcript") & (gtf.gene_biotype == "protein_coding")][cols_select].copy(deep=True)

    c = gtf1.groupby(['start','end']).seqname.count()
    # assert gtf1.shape[0] == np.sum(c) # sanity check
    print(f" duplicate promoters : {np.sum(c[c>1])} / {gtf1.shape[0]}")

    negstrand = gtf1.strand == "-"
    gtf1.loc[:,'start_n'] = gtf1.loc[:,'start']
    gtf1.loc[negstrand,'start_n'] = gtf1.loc[negstrand,'end']
    gtf1.start = gtf1.start_n - window//2
    gtf1.end = gtf1.start + window
    gtf1 = gtf1[cols_select]

    gtf1.start = gtf1.start.clip(lower=0)
    gtf1['name']=gtf1.apply(lambda rw:f"{rw.gene_id};{rw.gene_name};{rw.transcript_id};{rw.transcript_name};{rw.protein_id}",axis=1)
    gtf1.loc[pd.isna(gtf1.score),"score"] = "."
    gtf1[['seqname', 'start', 'end', 'name', 'score', 'strand']].to_csv(out_fname, sep = "\t", index=False, header=False)
    return out_fname

def process_promoter_v2(gtf_input, window, all_chroms, out_fname):
    """
    process gtf files for hg38 gencode 38, mm10 gencode vm 25
    """
    gtf = gtfparse.read_gtf(gtf_input)
    gtf = gtf[gtf['seqname'].apply(lambda x:x in all_chroms)]
    cols_select = ['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name', 'transcript_id', 'transcript_name', 'protein_id', 'score']
    gtf1 = gtf[(gtf.feature == "transcript") & (gtf.transcript_type == "protein_coding")][cols_select].copy(deep=True)

    c = gtf1.groupby(['start','end']).seqname.count()
    # assert gtf1.shape[0] == np.sum(c) # sanity check
    print(f" duplicate promoters : {np.sum(c[c>1])} / {gtf1.shape[0]}")

    negstrand = gtf1.strand == "-"
    gtf1.loc[:,'start_n'] = gtf1.loc[:,'start']
    gtf1.loc[negstrand,'start_n'] = gtf1.loc[negstrand,'end']
    gtf1.start = gtf1.start_n - window//2
    gtf1.end = gtf1.start + window
    gtf1 = gtf1[cols_select]

    gtf1.start = gtf1.start.clip(lower=0)
    gtf1['name']=gtf1.apply(lambda rw:f"{rw.gene_id};{rw.gene_name};{rw.transcript_id};{rw.transcript_name};{rw.protein_id}",axis=1)
    gtf1.loc[pd.isna(gtf1.score),"score"] = "."
    gtf1[['seqname', 'start', 'end', 'name', 'score', 'strand']].to_csv(out_fname, sep = "\t", index=False, header=False)
    return out_fname



def process_enhancer(enh_allfield, window, all_chroms, headers, out_path):
    # headers = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
    #        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
    #        'blockStarts']
    enh = pd.read_csv(enh_allfield, delimiter="\t", names=headers)
    enh = enh[enh.apply(lambda df:df['chrom'] in all_chroms, axis=1)]

    enh['chromStart'] = enh['thickStart']-window//2
    enh['chromStart'] = enh['chromStart'].clip(lower=0)
    enh['chromEnd'] = enh['chromStart'] + window

    enh.to_csv(out_path,sep='\t',header=False,index=False)
    return out_path

def splitCSV(input_path, out_paths, readArgs = {}, writeArgs = {}, prefix=None, split_num = None, suffix=None):
    """
    split an input csv file into N(=len(out_paths)) parts
    out_paths = [f"output_path_{i}.csv" for i in range(N)] #output path name function
    """
    N = len(out_paths)

    out_paths = []
    if N==0:
        for i in range(split_num):
            out_paths.append(prefix+str(i)+suffix)
        N = len(out_paths)

    df = pd.read_csv(input_path,**readArgs) # reading file
    max_lines = len(df)
    if max_lines==0:
        df.to_csv(out_paths[0],**writeArgs) 
        return 0
    numlines = int(np.floor(max_lines/N))
    low = np.arange(0,max_lines,numlines)[:N]
    high = np.concatenate((low[1:],[max_lines]))

    for i in range(N):
        df_new = df[low[i]:high[i]] # subsetting DataFrame based on index
        df_new.to_csv(out_paths[i],**writeArgs) # output file 
    return 0



def getChrDict(dna_out_file):
    """
    convert elements list into a dictionary keyed by chromosome; containing the (mid site, name ) tuple lists
    """
    def extractElem(elems):
        """for each element (promoters or enhancers) extract the chromosome info, start bp index, end bp index, and mid bp index"""
        matchObj = [re.match('chr([0-9XMY]+):([0-9]+)-([0-9]+)',pr) for pr in elems]
        Chr = [obj.group(1) for obj in matchObj]
        Start = [int(obj.group(2)) for obj in matchObj]
        End = [int(obj.group(3)) for obj in matchObj]
        Mid = [(x+y)//2 for x,y in zip(Start,End)]
        return Chr, Start, End, Mid

    elemData = np.load(dna_out_file,allow_pickle=True)
    elemChr, _, _, elemMid = extractElem(elemData['loc'])
    elemDF = pd.DataFrame({'chr':elemChr,'mid':elemMid, 'name':elemData['name']})
    # dictionary of (mid index, name) for array of elements ; dictionary keys = chromosomes; 
    # the arrays for each chromosome are sorted by mid
    return elemDF.groupby(['chr'])[['mid','name']].apply(lambda g: sorted(g.values.tolist(), key = lambda t: t[0])).to_dict()


def concat_PCHiC_PE(hicTSV,promoter_dna,enhancer_dna,selectCell='MK',threshold = 5, outputF=None, sampleFrac=None):
    """
    input : PCHiC tsv file, promoter .fa file, enhancer .fa file
    output : csv file with columns, baitPr, baitEnh, oePr, and oeEnh. 
        Corresponding to promoters and enhancers mid site intersecting with the bait and oe regions
    """
    print(selectCell)
    prChrDict =  getChrDict(promoter_dna)
    enhChrDict = getChrDict(enhancer_dna)

    pchicDF = pd.read_csv(hicTSV,delimiter="\t",dtype={'baitChr':str,'oeChr':str})
    pchicDF = pchicDF[pchicDF[selectCell]>=threshold]
    if sampleFrac:
        pchicDF = pchicDF.sample(frac=sampleFrac, axis=0)

    def intersectSortedElem(St,En,elemList):
        """
        returns all elements in elemList lying between St and En 
        i.e. element x in output array IFF St<=x<En
        """
        _elemList = np.array(elemList)[:,0].astype(np.int32)
        stIdx = bisect.bisect_left(_elemList,St)
        enIdx = bisect.bisect_left(_elemList,En)
        return elemList[stIdx:enIdx]

    def applyFunc(df,loci_type,elem_chr):
        Start = loci_type+'Start'
        End = loci_type+'End'
        chrom = loci_type+'Chr'
        # pdb.set_trace()
        if df[chrom] in elem_chr.keys():
            return intersectSortedElem(df[Start],df[End],elem_chr[df[chrom]])
        return []

    pchicDF['baitPr'] = pchicDF.apply(lambda df:applyFunc(df,'bait',prChrDict),axis=1) 
    pchicDF['baitEnh'] = pchicDF.apply(lambda df:applyFunc(df,'bait',enhChrDict),axis=1) 
    pchicDF['oePr'] = pchicDF.apply(lambda df:applyFunc(df,'oe',prChrDict),axis=1) 
    pchicDF['oeEnh'] = pchicDF.apply(lambda df:applyFunc(df,'oe',enhChrDict),axis=1) 

    if outputF:
        print(f"saving converted file to {outputF}..",end=" ",flush=True)
        pchicDF.to_pickle(outputF)
        # pd.read_pickle(outputF)
        print(".",flush=True)
    return pchicDF
