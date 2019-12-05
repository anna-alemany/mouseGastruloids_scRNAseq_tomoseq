# scRNA-seq scripts

### Summary of the data

We start with the merge dataset that contains 10X Genomics and all SORT-seq transcriptome data. Data is stored in mtx format for simplicity. 

All scRNA-seq datasets can be found in the Gene Expression Omnibus (GEO) under accession code [GSE123187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123187) and can be explored [here](https://avolab.hubrecht.eu/MouseGastruloids2019).

### Scripts repository subfolder

* **transcriptome_scanpy_10xSORT.py**

The [scanpy](https://scanpy.readthedocs.io/en/stable/) pipeline is used to perform the analysis of the scRNA-seq data. Data is filtered and ormalized, and variable genes are selected. Batch effects between different library-preparation technologies are corrected using [combat](https://www.ncbi.nlm.nih.gov/pubmed/16632515) and [Batch balanced KNN](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz625/5545955). 

All Leiden, Louvain and k-medoids clustering algorithms are tested, and differential gene expression analsysis for each approach is performed. 

* **cellcluster_annotation.py**

Based on the differential gene expression analysis results from the previous script, we select cell clusters that make sense with the biological properties of gastruloids. 

* **marioni_comparison.py**

Common genes between marker genes (P value < 0.01 and log2(fold-change) > 1.01) detected in the gastruloid cell clusters (identified with the previous script) and markers genes found for the different embryonic cell types defined in a [previously published mouse embryo scRNA-seq dataset](https://www.nature.com/articles/s41586-019-0933-9) are found in this script. P value for significance was assigned using a binomial test, where the probability of sharing a number of common marker genes between a gastruloid cell type and an embryonic cell type was determined by randomizing the list of marker genes for the embryonic cell type from the full list of marker genes in the embryonic cell types (n = 200). Comparisons to embryonic cell types found only at E8.5. and comparison to all embryonic cell types detected from E7.0 until E8.5. 

* **plotGeneUMAP.py**

This script provides a plot for the expression of any gene present in the dataset to the UMAP. 
