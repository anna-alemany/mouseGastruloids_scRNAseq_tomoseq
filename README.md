# Mouse gastruloids transcriptomics analysis

This repository contains the collection of scripts and python notebooks used to perform the transcriptomics analysis (scRNA-seq and tomo-seq) presented in [S. van den Brink et al](). 

The respository is organized in three subfolders:

#### 1. mapping_scripts

Contains all scripts necessary to map single-cell RNA sequencing data (scRNA-seq;, both 10X Genomics and SORT-seq approaches); as well as all the tomo-seq data (which is the same as SORT-seq in terms of mapping). 

#### 2. tomoseq_scripts

Contains all python notebooks required to filter and normalize tomo-seq data for mouse gastruloids (generated with two different cell lines), for E8.5 mouse embryos, and for the microarray microdisected PSM of E9.5 mouse embryos. Additionally, notebooks showing how the comparison between the different systems is done are also provided. 

#### 3. scRNAseq_scripts

Python scripts used to perform the analysis of the scRNA-seq data. Initially, scRNA-seq data from cells isolated from 5 days old gastruloids (cultured with two different cell lines) was generated using both 10X Genomics and SORT-seq protocols. Using the scanpy pipeline, all cells obtained with both technologies are analyzed togethers and batch effects are removed with combat and bbknn. Clustering and cell type calling is done using different algorithms and comparing our data to scRNA-seq data from cells isolated from E8.5 mouse embryos.  
