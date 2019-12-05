#!/usr/bin/env python3
import sys, os
from pandas.io.parsers import read_csv
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import sklearn.cluster as skcl
import scanpy as sc

option = '1'

sc.settings.figdir = '../results/manual_celltype_annotation/nodeCells/'

adata = sc.read_h5ad('../results/E14andLfng_v5/adata_clusters.h5ad')
df = read_csv('../results/manual_celltype_annotation/option1/umap_celltype_option1.tsv', sep = '\t', index_col=0)
adata.obs = adata.obs.merge(pd.DataFrame(df[['celltype']]).astype(str), how = 'inner', left_index = True, right_index = True)

tdf = pd.DataFrame(adata.raw.X.toarray())
tdf.index = adata.obs.index
tdf.columns = adata.raw.var.index
tdf = tdf.T

# differential gene expression of celltype groups
cluster = 'celltype'
sc.tl.rank_genes_groups(adata, cluster, method='t-test', n_genes = 100)
dex = {}
for cl in set(adata.obs[cluster]):
    keys = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
    dex[cl] = pd.DataFrame({k: adata.uns['rank_genes_groups'][k][cl] for k in keys})
    dex[cl] = dex[cl].set_index('names')
    dex[cl].to_csv('../results/manual_celltype_annotation/option1/dex_cl_' + cl.zfill(2) + '_option'+option+'.tsv', sep = '\t')

for cl in dex:
    genes = dex[cl].sort_values(by='logfoldchanges', ascending = False).index[:24]
    A = 3
    plt.figure(figsize = (A*1.6*4, A*6))
    for i, g in enumerate(genes):
        plt.subplot(6,4,i+1)
        cells = tdf.loc[g].sort_values().index
        plt.scatter(df.loc[cells,'u1'], df.loc[cells,'u2'], rasterized = True, s = 0.5, c = tdf.loc[g, cells])
        plt.colorbar(shrink = 0.5)
        plt.axis('off'); plt.title(g)
    plt.savefig(sc.settings.figdir + '/umap_cl'  + cl.zfill(2) + '_24topFC_' + cluster + '_option1.pdf')
    plt.close()
    
marker_genes = ['Hand2','Gata6','Meox1','Cxcl12','Bmp3','Meox2','Slc45a4','Uncx','Ripply2','Pcdh8',
        'Lfng','Tbx6','Hes7','Dll1','T','Fgf17','Wnt3a','Cdx2','Sox1','Pax6','Olfr959','Sox18','Igf1','Kdr',
        'Ahnak2','Tbx4','Hoxa11','Dppa5a','Nanog','Krt7','Epcam']

plt.close()

mfdf = pd.DataFrame({cl: {m: sum(tdf.loc[m, df[df['celltype']==cl].index]>0)/len(df[df['celltype']==cl]) for m in marker_genes} for cl in set(df['celltype'])})
mvdf = pd.DataFrame({cl: {m: tdf.loc[m, df[df['celltype']==cl].index].mean() for m in marker_genes} for cl in set(df['celltype'])})

plt.close()
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (4,6))
y = np.zeros(len(mvdf.columns))
x = range(len(mvdf.columns))
for m in marker_genes:
    im = ax.scatter(x, y, c = mvdf.loc[m], s = 50*mfdf.loc[m], cmap = 'bone_r', label = '_nolegend_')
    y += 1
py = 0.55*(len(marker_genes)+1)
px = 14
ax.scatter([px], [py], c = 'silver', s = 50, clip_on=False, label='100%')
ax.text(px+0.5, py-0.2, '100%')
ax.scatter([px], [py-1.5], c = 'silver', s = 25, clip_on=False, label='100%')
ax.text(px+0.5, py-1.5-0.2, '50%')
ax.scatter([px], [py-3], c = 'silver', s = 5, clip_on=False, label='100%')
ax.text(px+0.5, py-3-0.2, '10%')
ax.set_xticks(x); ax.set_xticklabels(mvdf.columns)
ax.set_yticks(range(len(marker_genes))); ax.set_yticklabels(marker_genes)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_xlabel('cluster')
cbaxes = fig.add_axes([0.9, 0.15,0.02, 0.3])
fig.colorbar(im, label = 'mean', shrink = 0.5, cax = cbaxes)
plt.savefig('../results/manual_celltype_annotation/option1/dotplot_markers.pdf', bbox_inches = 'tight')

# open marioni
import glob
files = glob.glob('/Users/anna/Dropbox/gastruloids/data/all-mGst-10xSCmergeSORTseq/results/manual_celltype_annotation/option1/dex_cl_*tsv')
dex = {f.rsplit('_')[4]: read_csv(f, sep = '\t', index_col=0) for f in files}

gene_df = read_csv('../data_mtx_format/genes.tsv', sep = '\t', index_col = 0, header = None, names = ['gene_name'])
#markers = read_csv('../../marioni/results/table_markerXcelltype_pooledTimePoints_aboveE7.tsv', sep = '\t', index_col=0)
markers = read_csv('../../marioni/results/table_markerXcelltype_E85.tsv', sep = '\t', index_col=0)
markers.index = [idx if idx not in gene_df.index else gene_df.loc[idx, 'gene_name'] for idx in markers.index]
markers = markers.loc[markers.index[markers.sum(axis=1)>0]]

# analysis
pvth = 0.1
for cl in dex:
    dx = dex[cl]
    dx = dx[(dx['pvals']<pvth)&(dx['logfoldchanges']>1.01)]
    dx = dx.sort_values(by='logfoldchanges', ascending=False)
    dx['rank'] = range(1,len(dx)+1)
    dex[cl] = dx

def findMatchingGenes(gene_names_markers, dex):
    rank = []
    for mg in gene_names_markers:
        if mg in dex.index:
            rank.append(dex.loc[mg,'rank'])
    return rank

cl2cell = pd.DataFrame(index = markers.columns, columns = dex.keys())
for k in markers.columns:
    for cl in dex.keys():
        r = findMatchingGenes(markers[k][markers[k]>0].index, dex[cl])
        cl2cell.loc[k, cl] = r

cl2cell_mean_df = cl2cell.applymap(lambda x: np.mean(x) if len(x) > 0 else 100)
cl2cell_len_df = cl2cell.applymap(lambda x: len(x))

#### binomial p-value ####
all_markers = markers.index

Nrep = 200
cl2cell_pv_df = pd.DataFrame(0, columns = dex.keys(), index = markers.columns)
for k in markers.columns:
    for cl in dex:
        l = []
        for nrep in range(Nrep):
            random_markers = np.random.choice(all_markers, size = markers[k].sum())
            l.append(len(findMatchingGenes(random_markers, dex[cl])))
        cnt = pd.Series(Counter(l))
        idxs = [idx for idx in cnt.index if idx >= len(cl2cell.loc[k, cl])]
        pv = cnt.loc[idxs].sum()/Nrep
        cl2cell_pv_df.loc[k,cl] = pv
        print(k, cl, pv)

cl2cell_pv_df.to_csv('../results/manual_celltype_annotation/option1/marioni_pvalues.tsv', sep = '\t')
cl2cell_len_df.to_csv('../results/manual_celltype_annotation/option1/marioni_numgenes.tsv', sep = '\t')
#cl2cell_pv_df.to_csv('../results/manual_celltype_annotation/option1/marioniE85_pvalues.tsv', sep = '\t')
#cl2cell_len_df.to_csv('../results/manual_celltype_annotation/option1/marioniE85_numgenes.tsv', sep = '\t')

cl2cell_pv_df = cl2cell_pv_df.loc[cl2cell_pv_df.index[(cl2cell_pv_df==1).sum(axis=1)!=len(cl2cell_pv_df.columns)]]
cl2cell = cl2cell.loc[cl2cell_pv_df.index]
cl2cell_len_df = cl2cell_len_df.loc[cl2cell_pv_df.index]

cl2cell_pv_fdf = cl2cell_pv_df # (cl2cell_pv_df <= 0.1)*cl2cell_pv_df + (cl2cell_pv_df > 0.1).astype(int)
cl2cell_pv_fdf = cl2cell_pv_df # cl2cell_pv_fdf.loc[cl2cell_pv_fdf.index[(cl2cell_pv_fdf==1).sum(axis=1) < len(cl2cell_pv_fdf.columns)]]

#### ranking scores ####
#th = 0.05
#th = 2
th = 1.# 0.2
rdfs = {}
for ct in cl2cell_pv_fdf.columns:
    print(ct)
    idxs = cl2cell_pv_fdf.loc[cl2cell_pv_fdf[cl2cell_pv_fdf[ct] < th].index,ct].sort_values().index
    ranks = cl2cell.loc[idxs, ct]
    ranks = ranks.apply(lambda x: sorted(x))
    probs = cl2cell_pv_fdf.loc[idxs,ct]
    gs = ranks.apply(lambda x: [dex[ct][dex[ct]['rank']==r].index[0] for r in x])
    rdf = pd.DataFrame({'pv': probs, 'ranks': ranks,'genes': gs})
    rdf['n_genes'] = [-len(rdf.loc[idx,'ranks']) for idx in rdf.index]
    rdf = rdf.sort_values(by=['pv','n_genes'])
    rdf['n_genes'] = -rdf['n_genes']
    rdfs[ct] = rdf
    rdfs[ct].to_csv('../results/manual_celltype_annotation/option1/marioniE85_topTargets_cl'+ct.zfill(2)+'.tsv', sep = '\t')
sort_cl = rdfs.keys()
sort_cl = ['1','2','3','4','5','6','7','8','9','10','11','12','13']

ms = list(set([ct for k in rdfs.keys() for ct in rdfs[k].index]))
ms_labels = {
        'Cardiomyocytes': 'Cardiomyocytes',
        'Mesenchyme': 'Mesenchyme',
        'Pharyngeal-mesoderm': 'Pharyngeal MS',
        'Allantois': 'Allantois',
        'Paraxial-mesoderm': 'Paraxial MS',
        'Notochord': 'Notochord',
        'Neural-crest': 'Neural crest',
        'Somitic-mesoderm': 'Somitic MS',
        'Caudal-Mesoderm': 'Caudal MS',
        'NMP': 'NMP',
        'Spinal-cord': 'Spinal cord',
        'Parietal-endoderm': 'Parietal endoderm',
        'Endothelium': 'Endothelium',
        'Haematoendothelial-progenitors': 'Haematoendo. prog.',
        'Blood-progenitors 1': 'Blood prog. 1',
        'Blood-progenitors 2': 'Blood-prog. 2',
        'ExE-endoderm': 'ExE endoderm',
        'PGC': 'PGC',
        'Forebrain-Midbrain-Hindbrain': 'Brain',
        'Visceral-endoderm': 'Visceral endoderm',
        'ExE-mesoderm': 'ExE mesoderm',
        'Gut': 'Gut', 'Rostral-neurectoderm': 'Rostral neurect.',
        'Surface-ectoderm': 'Surface ect.',
        'Intermediate-mesoderm': 'Intermediate MS',
        'Erythroid3': 'Erythroid3', 'Erythroid2': 'Erythroid2', 'Erythroid1': 'Erythroid1',
        'Def.-endoderm': 'Def. endoderm'
        }

ms_all = ['Cardiomyocytes','Mesenchyme','Pharyngeal-mesoderm','Allantois','Paraxial-mesoderm',
        'Notochord','Neural-crest', 'Somitic-mesoderm','Caudal-Mesoderm','NMP','Spinal-cord',
        'Parietal-endoderm', 'Endothelium','Haematoendothelial-progenitors', 'Blood-progenitors 1',
        'Blood-progenitors 2','ExE-endoderm','PGC','Forebrain-Midbrain-Hindbrain',
        'Visceral-endoderm',
        'ExE-mesoderm',
        'Gut', 'Rostral-neurectoderm',
        'Surface-ectoderm',
        'Intermediate-mesoderm', 'Erythroid3', 'Erythroid2', 'Erythroid1',
        'Def.-endoderm']
ms = [m for m in ms_all if m in ms]

plt.close()
A = 1.5
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (A*2.5,A*4)) # 3.8, 5.5
y = np.zeros(len(sort_cl))
x = range(len(sort_cl))
for i, m in enumerate(ms):
    im = ax.scatter(x, y, c = -np.log10(cl2cell_pv_fdf.loc[m, sort_cl]+1e-7), s =2*cl2cell_len_df.loc[m, sort_cl], cmap = 'Reds', label = '_nolegend_', vmax = 7)
    y += 1
y = len(ms)-0.4
ax.scatter([2], [y], c = 'silver', s = 60, clip_on=False, label='20')
ax.text(2.5, y-0.5, '30')
ax.scatter([6], [y], c = 'silver', s = 30, clip_on=False, label='10')
ax.text(6.5, y-0.5, '15')
ax.scatter([10], [y], c = 'silver', s = 10, clip_on=False, label='5')
ax.text(10.5, y-0.5, '5')
ax.set_xticks(x); ax.set_xticklabels(sort_cl)
ax.tick_params(axis = 'x', top = False, bottom = True, labelbottom = True, labeltop = False)
ax.set_yticks(range(len(ms))); ax.set_yticklabels([ms_labels[m] for m in ms])
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.spines['bottom'].set_visible(True)
ax.set_xlabel('cluster'); ax.xaxis.set_label_position('bottom')
cbaxes = fig.add_axes([-0.25, 0.11,0.3, 0.01])
fig.colorbar(im, label = '-log10(p-val)', shrink = 0.5, orientation = 'horizontal', cax = cbaxes)
plt.savefig('../results/manual_celltype_annotation/option1/dotplot_marioniComparison_pooled_option1_E85.pdf', bbox_inches = 'tight')


# fraction cell types in gastruloids
df['line'] = df.apply(lambda x: 'Lfng' if x['batch'] in ['10x','Lfng_sort'] else 'E14', axis = 1)

cntdf = pd.DataFrame({'Lfng': Counter(df[df['line']=='Lfng']['celltype']), 'E14': Counter(df[df['line']=='E14']['celltype'])})
fracdf = cntdf/cntdf.sum()

A = 1.2
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (A*3, A*1.6))
ax.bar(np.arange(len(fracdf.index))-0.2, fracdf.loc[fracdf.index,'E14'], label = 'E', width = 0.3, edgecolor='black')
ax.bar(np.arange(len(fracdf.index))+0.2, fracdf.loc[fracdf.index,'Lfng'], label = 'L', width = 0.3, edgecolor='black')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_xticks(x);
ax.set_xticklabels([])
#ax.set_xticklabels(fracdf.index)
ax.grid(which = 'major', axis = 'y', c = 'silver', lw = 0.3)
ax.set_axisbelow(True)
ax.set_xlim(-0.65,13.3)
ax.legend(frameon=False)
ax.set_ylabel('fraction');
#ax.set_xlabel('cluster')
plt.savefig(sc.settings.figdir + '/barplot_cellfraction_option1.pdf', bbox_inches ='tight')

# fraction cell types in embryo
mousedf = read_csv('../../marioni/SuppTable4_MetadataCellsAtlas.csv', index_col = 0)
mousedf = mousedf[mousedf['stage']=='E8.5']
cntdf = pd.Series(Counter(mousedf['celltype']))
cntdf.index = [str(c) for c in cntdf.index]
cntdf = cntdf.loc[[idx for idx in cntdf.index if idx != 'nan']]
fracdf = cntdf/cntdf.sum()
ms2 = [' '.join(m.rsplit('-')) for m in ms]
ms2 = [m if m != 'Forebrain Midbrain Hindbrain' else 'Forebrain/Midbrain/Hindbrain' for m in ms2]

plt.close()
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5.5, 1.4))
ax.bar(np.arange(len(ms2)), fracdf.loc[ms2[::-1]], width = 0.5, edgecolor='black', facecolor = 'silver')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_xticks(range(len(ms2)));
ax.set_xticklabels(['' for m in ms2[::-1]], rotation = 90)
ax.grid(which = 'major', axis = 'y', c = 'silver', lw = 0.3)
ax.set_axisbelow(True)
#ax.set_xlim(-1.6,len(ms)+1)
ax.set_xlim(-0.75, len(ms)-1+0.75)
ax.set_ylabel('fraction')
plt.savefig(sc.settings.figdir + '/barplot_cellfraction_mouseembryonic.pdf', bbox_inches ='tight')
