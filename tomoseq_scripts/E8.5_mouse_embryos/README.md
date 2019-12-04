# E8.5 mouse embryos

Two scripts are presented here:

* 1_Split-rawdata2mouseEmbryos-clean.ipynb

Here, barcodes that contain embryo material are selected from background for 3 replicates. Using expression of T/Bra, the posterior side of each embryo is identified.

* 2_SpatialPatterns_E85-batchCorrections-clean.ipynb

Here, we work with estimated transcripts count table (see mapping scripts to get more information). Spike-ins, ribosomal, and mitochondrial genes are removed for downstream analysis, together with Kcnq1ot1, Mir5109, Lars2, Malat1, Rn45s since these genes seem to be linked to mapping errors and have been shown to be erroneous in earlier studies.

Batch effects are corrected assuming by imposing the continuity of expression profiles along the anterior-posterior axis for each gene separately. 

Next (in the same script), data in each embryo is normalized to the median number of unique transcripts per slide and the z-score profile along the sections for each gene is extracted.

Gene reproducibility analysis between replicates is performed by computing the Pearson correlation coefficient between the AP expression pattern (in z-score units) of two different samples for all possible pairs of replicates. Linearly interpolated gene expression profiles are used when the number of sections is different between replicates. To assess for significant correlations, we randomly generate 5,000 expression profiles with the same number of sections as in the pair of replicates and determine a threshold for the correlation value at which less than n random profiles have larger correlation values (n = 100 for P value < 0.01; n = 500 for P value < 0.05, etc; Supplementary Table 5). Adjusted P values are obtained with the Benjamini/Hochberg correction.

Only genes that are significantly correlated (P value < 0.01) in at least five possible pairs of replicates are considered as reproducible between replicates. In the script, reproducible genes are clustered using self-organizing maps and hierarchical clsutering according to their expression profile.

© 2019 
