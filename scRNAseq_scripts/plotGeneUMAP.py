#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
import matplotlib.pyplot as plt
import scanpy as sc

try:
    g = sys.argv[1]
except:
    sys.exit('Please, give input gene name')

sc.settings.figdir = '../results/manual_celltype_annotation/option1/umap_geneExpression/'
adata = sc.read_h5ad('../results/E14andLfng_v5/adata_clusters.h5ad')
df = read_csv('../results/manual_celltype_annotation/option1/umap_celltype_option1.tsv', sep = '\t', index_col=0)

tdf = pd.DataFrame(adata.raw.X.toarray())
tdf.index = adata.obs.index
tdf.columns = adata.raw.var.index
tdf = tdf.T

fig, ax = plt.subplots(ncols = 1, nrows = 1, figsize = (5*1.5, 5))
cells = tdf.loc[g].sort_values(ascending=True).index
im = ax.scatter(df.loc[cells,'u1'], df.loc[cells,'u2'], c = tdf.loc[g, cells], s = 0.5, cmap = 'RdYlBu_r', rasterized = True)
ax.set_xticks([]); ax.set_yticks([])
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False); ax.spines['left'].set_visible(False)
ax.set_title(g)
plt.savefig(sc.settings.figdir + '/umap-' + g + '.pdf', bbox_inches = 'tight')
