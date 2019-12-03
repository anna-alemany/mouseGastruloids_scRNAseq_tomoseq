# Mapping scripts

## Introduction

### scRNA-seq (SORT-seq and 10x Genomics). 

For scRNA-seq, cells extracted from 120 h gastruloids were processed using either SORT-seq 
(CEL-seq2 based scRNA-seq on cells that were sorted into 384-well plates32) or using the 10x Genomics Chromium Single Cell 3' (v3 Chemistry)
gene expression kit, according to manufacturer's instructions. 

### Tomo-seq

Tomo-seq was performed using a robotized (SORT-seq-based) version of a previously published tomo-seq protocol. 
Briefly, 120 h gastruloids or E8.5 mouse embryos were embedded in cryosolution (Leica, 14020108926), snap-frozen on dry-ice, 
stored at -80 °C and sectioned using a cryotome. 

Sections were collected in the wells of a Hard-Shell PCR Low-profile, semi-skirted 96-well plate (Bio-rad, HSL9601) that was already 
prefilled with mineral oil (Sigma, M8410-1L) and CEL-seq2 primers. 
For each well, a unique, barcoded CEL-seq2 primer was used, which allowed us to pool the content of the wells after second strand synthesis.
To sequence the mRNA content of the wells, SORT-seq (robotized CEL-seq2 based scRNA-seq) was performed using a Nanodrop II liquid 
handling platform (GC biotech).

### Sequencing

Sequencing was performed on the Illumina Next-seq sequencing platform. 
For SORT-seq and tomo-seq, paired end (75 bp) sequencing was performed; 
for 10x Genomics, sequencing was performed according to 10x Genomics manufacturer’s instructions (Read1, 28 cycles; Index i7, 8 cycles; Read2, 91 cycles).

### Mapping sequencing data. 

For SORT-seq and tomo-seq, the first 6 bases of read 1 contain the unique molecular identifier (UMI) and the next 7 bases contain the cell or section barcode. 

For 10x Genomics, the first 16 bases of read 1 contain the cell barcode, and the next 12 contain the UMI. 

For all sequencing experiments, read 2 contains the biological information. 
Reads 2 with a valid cell/section barcode were selected, 
trimmed using TrimGalore (v0.4.3) with default parameters, and mapped using STAR (v2.5.3a) 
with default parameters to the mouse mm10 genome (Ensembl 93). 
Only reads mapping to gene bodies (exons or introns) were used for downstream analysis. 

Reads mapping simultaneously to an exon and to an intron were assigned to the exon. 
For each cell or section, the number of transcripts was obtained as previously described (D.Grun et al). 
We refer to transcripts as unique molecules based on UMI correction. 

## Description of scripts

### SORT-seq and tomo-seq

Let's assume we start with the following fastq files:

| Read 1 | Read 2 |
| --- | --- |
| library_L001_R1_001.fastq.gz | library_L001_R2_001.fastq.gz |
| library_L002_R1_001.fastq.gz | library_L002_R2_001.fastq.gz |
| library_L003_R1_001.fastq.gz | library_L003_R2_001.fastq.gz | 
| library_L004_R1_001.fastq.gz | library_L004_R2_001.fastq.gz | 

To map them to the mouse genome, we need to type in the terminal:
```{bash}
submit_array_starmap.sh library_L00 output_name
```
The script *submit_array_starmap.sh* will submit a sequence of scripts that perform the required steps to map our library to the mouse genome and produce the count tables.

In order for the script to run, we need TrimGalore-0.4.3, cutadapt, STAR, python3, together with all the scripts called.

This will produce a total of 9 files:

| | spliced | unspliced | total | 
| --- | --- | --- | --- | 
**read counts** | output_name_spliced.coutc.tsv | output_name_unspliced.coutc.tsv | output_name_total.coutc.tsv | 
**observed UMIs** | output_name_spliced.coutb.tsv | output_name_spliced.coutb.tsv | output_name_total.coutb.tsv | 
**estimated transcripts** | output_name_unspliced.coutt.tsv | output_name_spliced.coutt.tsv | output_name_total.coutt.tsv | 

Unspliced, spliced or total denotes whether the read contains some region in an intron (unspliced) or an exon (spliced) of the annotated gene. Total does not take introns/exons into account. coutc refers to the total number of reads, coutb to the total number of observed UMIs, and coutt is the total number of unique transcripts. The last one is obtained from coutb by applying the Poisson correction described by D. Grun in his Nature paper.

