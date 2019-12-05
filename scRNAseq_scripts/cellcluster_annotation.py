#!/usr/bin/env python3
import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing

mycolors = ["#008941", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#FFB500", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FF0000", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09"]

umap = read_csv('../results/E14andLfng_v5/umap.tsv', sep = '\t', index_col=0)
umap.columns = ['u1','u2']
df = read_csv('../results/E14andLfng_v5/obs_info.tsv', sep = '\t', index_col=0)

df = umap.merge(df, how = 'inner', left_index = True, right_index = True)
df.to_csv('../results/manual_celltype_annotation/umap-clusters.tsv', sep = '\t')

#### cluster plots
cs = {b: i+1 for i, b in enumerate(set(df['batch']))}
cells = list(df.index)
while len(set([df.loc[c,'batch'] for c in cells[:3]])) != 3:
    np.random.shuffle(cells)
plt.scatter(df.loc[cells[0],'u1'], df.loc[cells[0],'u2'], s = 2, color = mycolors[cs[df.loc[cells[0],'batch']]], rasterized = True, alpha = 0.5, label = df.loc[cells[0],'batch'])
plt.scatter(df.loc[cells[1],'u1'], df.loc[cells[1],'u2'], s = 2, color = mycolors[cs[df.loc[cells[1],'batch']]], rasterized = True, alpha = 0.5, label = df.loc[cells[1],'batch'])
plt.scatter(df.loc[cells[2],'u1'], df.loc[cells[2],'u2'], s = 2, color = mycolors[cs[df.loc[cells[2],'batch']]], rasterized = True, alpha = 0.5, label = df.loc[cells[2],'batch'])
plt.scatter(df.loc[cells,'u1'], df.loc[cells,'u2'], s = 2, color = [mycolors[cs[df.loc[idx,'batch']]] for idx in cells], rasterized = True, alpha = 0.5, label = '_nolabel_')
plt.legend(loc = 2, bbox_to_anchor = (1,1))
plt.xlabel('umap 1'); plt.ylabel('umap 2')
plt.title('Batch')
plt.savefig('../results/manual_celltype_annotation/umap_batch.pdf', bbox_inches = 'tight')

for i, cl in enumerate(set(df['leiden'])):
    cells = df[df['leiden']==cl].index
    plt.scatter(df.loc[cells,'u1'], df.loc[cells,'u2'], s = 2, color = mycolors[i], label = cl, rasterized = True)
plt.legend(loc = 2, bbox_to_anchor = (1,1))
plt.xlabel('umap 1'); plt.ylabel('umap 2')
plt.title('Leiden clustering')
plt.savefig('../results/manual_celltype_annotation/umap_leiden.pdf', bbox_inches = 'tight')
plt.close()

for i, cl in enumerate(set(df['kmedoids11'])):
    cells = df[df['kmedoids11']==cl].index
    plt.scatter(df.loc[cells,'u1'], df.loc[cells,'u2'], s = 2, color = mycolors[i], label = cl, rasterized = True)
plt.legend(loc = 2, bbox_to_anchor = (1,1))
plt.xlabel('umap 1'); plt.ylabel('umap 2')
plt.title('k-medoids clustering')
plt.savefig('../results/manual_celltype_annotation/umap_kmedoids11.pdf', bbox_inches = 'tight')
plt.close()

#### combination of leiden and k-medoids
## option 1
df['celltype'] = 'nan'
for idx in df.index:
    if df.loc[idx, 'leiden'] == 8 and df.loc[idx, 'u2'] > 2.5 and df.loc[idx, 'celltype'] == 'nan':
        df.loc[idx, 'celltype'] = 'L08'
    if df.loc[idx, 'leiden'] == 10 and df.loc[idx, 'celltype'] == 'nan':
        df.loc[idx, 'celltype'] = 'L10'
    if df.loc[idx, 'leiden'] == 9 and df.loc[idx, 'celltype'] == 'nan':
        df.loc[idx, 'celltype'] = 'L09'
    if df.loc[idx, 'leiden'] == 11 and df.loc[idx, 'celltype'] == 'nan':
        df.loc[idx, 'celltype'] = 'L11'
    if df.loc[idx, 'kmedoids11'] == 3 and df.loc[idx, 'celltype'] == 'nan' and df.loc[idx, 'u2'] > 2.5:
        df.loc[idx, 'celltype'] = 'K03'
for idx in df.index:
    if df.loc[idx, 'celltype'] == 'nan':
        df.loc[idx, 'celltype'] = 'LL' + str(df.loc[idx, 'leiden']).zfill(2)
    if df.loc[idx, 'celltype'] == 'LL12':
        df.loc[idx, 'celltype'] = 'LL00'

rename = {'K03': 1, 'LL07': 2, 'LL04': 3, 'LL00': 3, 'LL02': 4,
          'LL06': 5, 'LL05': 6, 'LL03': 7, 'LL01': 8, 'LL08': 9,
          'L08': 10, 'L10': 11, 'L09': 12, 'L11': 13}

df['celltype'] = df.apply(lambda x: rename[x['celltype']], axis = 1)

df.to_csv('../results/manual_celltype_annotation/option1/umap_celltype_option.tsv', sep = '\t')

plt.close()
A = 4
plt.figure(figsize = (A*1.3,A))
ax = plt.subplot(111)
for i, cl in enumerate(set(df['celltype'])):
    cells = df[df['celltype']==cl].index
    ax.scatter(df.loc[cells,'u1'], df.loc[cells,'u2'], s = 1.5, color = mycolors[i], label = cl, rasterized = True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.legend(loc = 3, ncol = 13)#, bbox_to_anchor = (1,1))
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol = 5)
ax.text(-6,9.5,'n=25202')
ax.set_xlabel('umap 1'); ax.set_ylabel('umap 2')
ax.set_xticks([]); ax.set_yticks([])
plt.savefig('../results/manual_celltype_annotation/option1/umap_celltype_option.pdf', bbox_inches = 'tight')
plt.close()
