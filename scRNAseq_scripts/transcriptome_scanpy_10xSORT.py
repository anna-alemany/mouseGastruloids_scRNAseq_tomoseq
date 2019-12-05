#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
import scanpy as sc
import bbknn

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.figdir = '../results/E14andLfng_v5_posREVISION2/'
sc.settings.set_figure_params(dpi=80)

#### Read data for scanpy ####
adata = sc.read_10x_mtx(
        '../data_mtx_format',
        var_names = 'gene_symbols',
        cache = True)

#### Pre-processing ####
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
adata.obs['n_genes'] = (adata.X>0).sum(axis=1).A1
adata.obs['protocol'] = [idx.rsplit('_')[-1] for idx in adata.obs.index]
adata.obs['experiment'] = ['_'.join(idx.rsplit('_')[-2:]) for idx in adata.obs.index]
adata.obs['experiment'] = [es if 'sort' in es else es.rsplit('-')[-1] for es in adata.obs['experiment']]
adata.obs['batch'] = [es if '10x' not in es else '10x' for es in adata.obs['experiment']]
adata.obs['plate'] = ['-'.join(idx.rsplit('-')[1:]) for idx in adata.obs.index]

sc.pl.scatter(adata, x = 'n_counts', y = 'n_genes', color = 'experiment', save = '_genesVScounts_experiment.pdf')
sc.pl.scatter(adata, x = 'n_counts', y = 'n_genes', color = 'protocol', save = '_genesVScounts_protocol.pdf')
sc.pl.scatter(adata, x = 'n_counts', y = 'n_genes', color = 'batch', save = '_genesVScounts_batch.pdf')

sc.pp.filter_genes(adata, min_cells = 3) # filtered out 2030 genes that are detected in less than 3 cells
adata = adata[adata.obs['n_genes'] < 8000, :]
adata = adata[adata.obs['n_counts'] < 40000, :]
sc.pl.scatter(adata, x = 'n_counts', y = 'n_genes', color = 'experiment', save = '_genesVScounts_experiment_postfilter.pdf')

adata = adata[adata.obs['n_counts'] > 1000, :]
adata = adata[adata.obs['n_genes'] > 700, :]

sc.pl.scatter(adata, x = 'n_counts', y = 'n_genes', color = 'experiment', save = '_genesVScounts_experiment_postfilter.pdf')
print(adata.shape)

#### normalization ####
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata) # X = \log(X + 1)`, where :math:`log` denotes the natural logarithm.
adata.raw = adata

#### gene variability ####
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=5, min_disp=0.5, n_bins = 30)
sc.pl.highly_variable_genes(adata, save = True)

adata = adata[:, adata.var['highly_variable']]

#### batch corrections ####
sc.pp.combat(adata, key='batch')

#### PCA ####
sc.pp.regress_out(adata, ['n_counts'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps = 100)

sc.pl.pca(adata, color=['experiment'], save = '_10XvsSort.pdf')
sc.pl.pca(adata, color=['Meox1','T','Mt1', 'Gata6', 'Tmsb4x', 'Crabp1', 'Sox17', 'Pax6'], save = '_genes.pdf')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50, save = True)

#### neighborhood graph ####
npca = 40
#sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = npca, random_state = 1971723, metric = 'manhattan')
bbknn.bbknn(adata, batch_key = 'batch', n_pcs = npca, metric = 'manhattan')

#### UMAP or tsne ####
sc.tl.umap(adata, n_components = 2, random_state = 924163)#, spread = 1, min_dist = 0.75)
sc.pl.umap(adata, color = ['batch'], save = '_batch.pdf')
sc.pl.umap(adata, color=['Meox1','T','Mt1', 'Gata6', 'Tmsb4x', 'Cers1', 'Sox17', 'Pax6'], save = '_genes.pdf')

sc.tl.tsne(adata, random_state = 924163)#, spread = 1, min_dist = 0.75)
sc.pl.tsne(adata, color = ['batch'], save = '_batch.pdf')
sc.pl.tsne(adata, color=['Meox1','T','Mt1', 'Gata6', 'Tmsb4x', 'Cers1', 'Sox17', 'Pax6'], save = '_genes.pdf')

#### Clustering ####
sc.tl.leiden(adata, random_state = 31029)
sc.tl.louvain(adata, random_state = 3381)

df = pd.DataFrame(adata.obsm['X_pca'], index = adata.obs.index)
dist_df = 1 - df[range(npca)].T.corr()

def kmedoids(dist, numclusters):
    func_min = lambda x, medoids: medoids[dist.loc[medoids,x] == dist.loc[medoids,x].min()]
    findMedoid = lambda idxs, dists: pd.Series({idx: dist.loc[idx, idxs].sum() for idx in idxs}).sort_values().index[0]

    vj = pd.Series()
    for col in dist.columns:
        c = dist.loc[col].sum()
        vj.loc[col] = dist[col].sum()/c

    medoids = np.array(vj.sort_values(ascending=True).index[:numclusters])
    clusters = pd.Series({col: func_min(col, medoids)[0] for col in dist.columns})
    D0 = sum([dist.loc[col, clusters[col]] for col in clusters.index])

    while 2>1:
        medoids_new = []
        for m in medoids:
            idxs = clusters[clusters == m].index
            medoids_new.append(findMedoid(idxs, dist))
        medoids_new = np.array(medoids_new)

        clusters_new = pd.Series({col: func_min(col, medoids_new)[0] for col in dist.columns})

        D1 = sum([dist.loc[col, clusters_new[col]] for col in clusters_new.index])

        medoids = medoids_new
        clusters = clusters_new
        print(D0,D1)

        if D1 < D0:
            D0 = D1
        else:
            break
    return medoids, clusters

medoids, kmedoids_clusters = kmedoids(dist_df, numclusters = 9)
medoids = {m: i for i, m in enumerate(medoids)}
adata.uns['kmedoids09'] = {'params': {'numclusters' : 9}}
adata.obs['kmedoids09'] = [str(medoids[kmedoids_clusters[idx]]) for idx in adata.obs.index]

medoids, kmedoids_clusters = kmedoids(dist_df, numclusters = 10)
medoids = {m: i for i, m in enumerate(medoids)}
adata.uns['kmedoids10'] = {'params': {'numclusters' : 10}}
adata.obs['kmedoids10'] = [str(medoids[kmedoids_clusters[idx]]) for idx in adata.obs.index]

medoids, kmedoids_clusters = kmedoids(dist_df, numclusters = 11)
medoids = {m: i for i, m in enumerate(medoids)}
adata.uns['kmedoids11'] = {'params': {'numclusters' : 11}}
adata.obs['kmedoids11'] = [str(medoids[kmedoids_clusters[idx]]) for idx in adata.obs.index]

medoids, kmedoids_clusters = kmedoids(dist_df, numclusters = 12)
medoids = {m: i for i, m in enumerate(medoids)}
adata.uns['kmedoids12'] = {'params': {'numclusters' : 12}}
adata.obs['kmedoids12'] = [str(medoids[kmedoids_clusters[idx]]) for idx in adata.obs.index]

medoids, kmedoids_clusters = kmedoids(dist_df, numclusters = 13)
medoids = {m: i for i, m in enumerate(medoids)}
adata.uns['kmedoids13'] = {'params': {'numclusters' : 13}}
adata.obs['kmedoids13'] = [str(medoids[kmedoids_clusters[idx]]) for idx in adata.obs.index]

medoids, kmedoids_clusters = kmedoids(dist_df, numclusters = 14)
medoids = {m: i for i, m in enumerate(medoids)}
adata.uns['kmedoids14'] = {'params': {'numclusters' : 14}}
adata.obs['kmedoids14'] = [str(medoids[kmedoids_clusters[idx]]) for idx in adata.obs.index]

sc.pl.umap(adata, color = ['leiden', 'louvain', 'kmedoids09'], save = '_clusters09.pdf')
sc.pl.umap(adata, color = ['leiden', 'louvain', 'kmedoids10'], save = '_clusters10.pdf')
sc.pl.umap(adata, color = ['leiden', 'louvain', 'kmedoids11'], save = '_clusters11.pdf')
sc.pl.umap(adata, color = ['leiden', 'louvain', 'kmedoids12'], save = '_clusters12.pdf')
sc.pl.umap(adata, color = ['leiden', 'louvain', 'kmedoids13'], save = '_clusters13.pdf')
sc.pl.umap(adata, color = ['leiden', 'louvain', 'kmedoids14'], save = '_clusters14.pdf')

sc.pl.umap(adata, color = ['kmedoids09','kmedoids10','kmedoids11','kmedoids12','kmedoids13','kmedoids14'], save = '_kmedoids.pdf')

adata.obs.to_csv('../results/E14andLfng_v5/obs_info.tsv', sep = '\t')
adata.var.to_csv('../results/E14andLfng_v5/obs_var.tsv', sep = '\t')

umap_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index)
umap_df.to_csv(sc.settings.figdir + 'umap.tsv', sep = '\t')
tsne_df = pd.DataFrame(adata.obsm['X_tsne'], index=adata.obs.index)
tsne_df.to_csv(sc.settings.figdir + 'tsne.tsv', sep = '\t')

adata.write(sc.settings.figdir + 'adata_clusters.h5ad', compression='gzip')

#### differential gene expresssion analysis ####
dex = {}
for cluster in ['leiden','louvain', 'kmedoids11']:
    sc.tl.rank_genes_groups(adata, cluster, method='t-test')
    dex[cluster] = {}
    for cl in set(adata.obs[cluster]):
        keys = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
        dex[cluster][cl] = pd.DataFrame({k: adata.uns['rank_genes_groups'][k][cl] for k in keys})
        dex[cluster][cl] = dex[cluster][cl].set_index('names')
        dex[cluster][cl].to_csv('../results/E14andLfng_v5/dex_cl' + cl.zfill(2) + '_' + cluster + '.tsv', sep = '\t')

    for cl in dex[cluster]:
        sc.pl.umap(adata, color = dex[cluster][cl].sort_values(by='logfoldchanges', ascending = False).index[:24], save = '_cl' + cl.zfill(2) + '_24topFC_' + cluster + '.pdf')
        
#### Interesting genes ###
# somitogenesis wavefront
sc.pl.umap(adata, color = ['Mesp1', 'Mesp2', 'Ripply2', 'Lfng', 'Notch1'], save = '_wavefront.pdf')
# cranial mesoderm
sc.pl.umap(adata, color = ['Mest', 'Cxcl12', 'Cdh11', 'Eya2'], save = '_cmesoderm.pdf')
# somites
sc.pl.umap(adata, color = ['Tcf15','Aldh1a2','Uncx','Tbx18'], save = '_somiticMesoderm.pdf')
# neural tissue
sc.pl.umap(adata, color = ['Sox1', 'Sox2','Foxg1','Pax6'], save = '_neural.pdf')
# tailbud
sc.pl.umap(adata, color = ['T','Cdx2','Wnt3a','Fgf8'], save = '_PSM.pdf')
# head mesenchyme / neural crest
sc.pl.umap(adata, color = ['Crabp1','Btbd17','Elavl3','Elavl4','Robo1','Celsr3','Onecut2'], save = '_neuralCrest.pdf')
# heart
sc.pl.umap(adata, color = ['Hand1','Hand2','Gata4','Gata6','Bmp4'], save = '_heart.pdf')
# allantois
sc.pl.umap(adata, color = ['Tbx4','Hoxa11','Krt8','Krt18','Bmp7'], save = '_allantois.pdf')
# node (pit)
sc.pl.umap(adata, color = ['Foxa1','Foxa2','Sox17','Epcam','Prdm1','Trh'], save = '_nodePit.pdf')
# node (crown)
sc.pl.umap(adata, color = ['T','Nodal','Lefty1','Lefty2','Cdh1','Trh'], save = '_nodeCrown.pdf')
# late endoderm
sc.pl.umap(adata, color = ['Mixl1','Eomes','Gsc','Lefty1','Lefty2'], save = '_lateEndoderm.pdf')
# early endoderm
sc.pl.umap(adata, color = ['Prdm1','Irx1'], save = '_earlyEndoderm.pdf')
# endothelial
sc.pl.umap(adata, color = ['Kdr','Sox7','Sox17','Tie1','Cdh5'], save = '_endothelial.pdf')
# PGCs
sc.pl.umap(adata, color = ['Mt1', 'Mt2', 'Dppa5a', 'Esrrb', 'Otx2'], save = '_PGC.pdf')
