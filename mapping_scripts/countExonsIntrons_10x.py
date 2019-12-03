#!/usr/bin/env python3
import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
from collections import Counter

try:
    bedintronfile=sys.argv[1]
    bedexonfile = sys.argv[2]
    cellBarcodes = sys.argv[3]
    output = sys.argv[4]
except:
    sys.exit("Please, provide input bed (1) intron and (2) exon files; (3) cell barcodes; (4) output file")

cells = pd.read_csv(cellBarcodes, header = None)
cells = list(cells[0])

cnt = {}
cntreadi = {}
cntreade = {}
cntread = {}
for label, bedfile in [('intron', bedintronfile), ('exon', bedexonfile)]:
    with open(bedfile) as f:
        for line in f:
            ch, x0, x1, name, strand, gene = line.rsplit('\t')
            gene = '__'.join([gene.rstrip(), ch])
            umi = name.rsplit(':')[8]
            cell = name.rsplit(':')[7]
            if cell in cells:
                try:
                    cnt[cell][gene][umi].update([label])
                except:
                    try:
                        cnt[cell][gene][umi] =  Counter([label])
                    except:
                        try:
                            cnt[cell][gene] = {umi: Counter([label])}
                        except:
                            cnt[cell] = {gene: {umi: Counter([label])}}
                if label == 'intron':
                    try:
                        cntreadi[cell][gene] += 1
                    except:
                        try:
                            cntreadi[cell][gene] =  1
                        except:
                            cntreadi[cell] = {gene: 1}

                if label == 'exon':
                    try:
                        cntreade[cell][gene] += 1
                    except:
                        try:
                            cntreade[cell][gene] = 1
                        except:
                            cntreade[cell] = {gene: 1}

                try:
                    cntread[cell][gene] += 1
                except:
                    try:
                        cntread[cell][gene] = 1
                    except:
                        cntread[cell] = {gene: 1}

df = pd.DataFrame(cnt)
dfri = pd.DataFrame(cntreadi)
dfre = pd.DataFrame(cntreade)
dfr = pd.DataFrame(cntread)

cols = sorted(df.columns)
df = df[cols]
cols = sorted(dfri.columns)
dfri = dfri[cols].fillna(0).astype(int)
cols = sorted(dfre.columns)
dfre = dfre[cols].fillna(0).astype(int)
cols = sorted(dfr.columns)
dfr = dfr[cols].fillna(0).astype(int)

dfri.to_csv(output + '_unspliced.coutc.tsv', sep = '\t')
dfre.to_csv(output + '_spliced.coutc.tsv', sep = '\t')
dfr.to_csv(output + '_total.coutc.tsv', sep = '\t')

def countUnsplicedMolecules(x):
    y = 0
    if type(x) == dict:
        for umi in x:
            if 'intron' in x[umi]:
                y += 1
    return y

def countSplicedMolecules(x):
    y = 0
    if type(x) == dict:
        for umi in x:
            if 'intron' not in x[umi]:
                y += 1
    return y

def bc2trans(x):
    if x >= K:
        t = np.log(1.-(float(K)-1e-3)/K)/np.log(1.-1./K)
    elif x > 0 and x < K:
        t = np.log(1.-float(x)/K)/np.log(1.-1./K)
    elif x == 0:
        t = 0
    return  t

udf = df.applymap(countUnsplicedMolecules)
sdf = df.applymap(countSplicedMolecules)
bdf = df.applymap(lambda x: len(x) if type(x)==dict else 0)

udf.to_csv(output + '_unspliced.coutb.tsv', sep = '\t')
sdf.to_csv(output + '_spliced.coutb.tsv', sep = '\t')
bdf.to_csv(output + '_total.coutb.tsv', sep = '\t')

K = 4**len(umi)

tudf = udf.applymap(bc2trans)
tsdf = sdf.applymap(bc2trans)
tbdf = bdf.applymap(bc2trans)

tudf.to_csv(output + '_unspliced.coutt.tsv', sep = '\t')
tsdf.to_csv(output + '_spliced.coutt.tsv', sep = '\t')
tbdf.to_csv(output + '_total.coutt.tsv', sep = '\t')
