# E14 mouse gastruloids

Two python notebooks are presented here: 

* **1_Split-rawdata2gastruloids-clean.ipynb**

Here, barcodes that contain gastruloid material are selected for each of the 5 replicates. Using expression of T/Bra, the posterior side of each gastruloid is identified. 

* **2_similarities_5dAA_clean.ipynb**

Here, we work with estimated transcripts count table (see mapping scripts to get more information). Spike-ins, ribosomal, and mitochondrial genes are removed for downstream analysis, together with _Kcnq1ot1_, _Mir5109_, _Lars2_, _Malat1_, _Rn45s_ since these genes seem to be linked to mapping errors and have been shown to be erroneous in earlier studies.

Next (in the same script), data in each gastruloid is normalized to the median number of unique transcripts per slide and the z-score profile along the sections for each gene is extracted. 

Gene reproducibility analysis between replicates is performed by computing the Pearson correlation coefficient between the AP expression pattern (in z-score units) of two different samples for all possible pairs of replicates. Linearly interpolated gene expression profiles are used when the number of sections is different between replicates. To assess for significant correlations, we randomly generate 5,000 expression profiles with the same number of sections as in the pair of replicates and determine a threshold for the correlation value at which less than n random profiles have larger correlation values (n = 100 for P value < 0.01; n = 500 for P value < 0.05, etc; Supplementary Table 5). Adjusted P values are obtained with the Benjamini/Hochberg correction. 

Only genes that are significantly correlated (P value < 0.01) in at least five possible pairs of replicates are considered as reproducible between replicates. In the script, reproducible genes are clustered using self-organizing maps and hierarchical clustering according to their expression profile. 
