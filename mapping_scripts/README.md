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

Let's assume we start with the following fastq files (common scheme in both 10X and SORT-seq and tomo-seq):

| Read 1 | Read 2 |
| --- | --- |
| library_L001_R1_001.fastq.gz | library_L001_R2_001.fastq.gz |
| library_L002_R1_001.fastq.gz | library_L002_R2_001.fastq.gz |
| library_L003_R1_001.fastq.gz | library_L003_R2_001.fastq.gz | 
| library_L004_R1_001.fastq.gz | library_L004_R2_001.fastq.gz | 

### SORT-seq and tomo-seq data

Here, the number of cells per library is well known since transcripts from each cells have been barcoded via IVT. 

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

### 10X Genomics data

(because we wanted to keep the mapping as similar as possible to SORTseq, we did not use cellranger and custom-made code was written).

In this case, the number of cells is not know. To identify the number of cells detected in our 10X library, we first run the script: 
```{bash}
submit_array_starmap10X.sh library_L00 output_name
```
which, as before, will merge the lanes, extract cell barcodes and UMIs from each read, trim the reads and map. Before generating any count table, the script will count the number of reads per barcode. A file called **output_name_cbc_trimmed_rawCELLS.txt** will be generated, which reads as: 
```{bash}
  16943 TATCGCCCACAGCTGC
  16593 GTACAACGTAGACGTG
  16532 TCTGGCTCACACCGCA
  16080 CATCCCAGTAGCGCTC
  15994 AAAGGTAAGCCTGAGA
  15502 TCATGCCTCACTCCGT
  15124 GCAGCCAAGACCATGG
  14928 TCCATCGTCTCCTGCA
  14878 GGATGTTTCGGCCTTT
  14844 GACCAATCAACCACAT
  etc
```
where the first column indicates the number of reads and the second the number of barcodes (both depend on the library). Exploring the number of reads per barcode (see figure below) it becomes natural to identify approximately the number of cells that have been sequenced. In the example shown here, the number of cells (vertical line) is set to 4,768.

![alt text](https://github.com/anna-alemany/mouseGastruloids_scRNAseq_tomoseq/blob/master/mapping_scripts/cells.jpg "Reads per Barcode")

Now, we create a text file containing the barcodes that are assigned to valid cells. For example, typing:
```{bash}
head -4768 output_name_cbc_trimmed_rawCELLS.txt | awk '{print $2}' > output_name_cbc_trimmed_selectedCELLS.txt 
```

Finally, to produce the desired count table with transcriptome information, we use:
```{bash}
./submit_coutExonsIntrons10x.sh output_name_cbc_trimmed_starAligned.sortedByCoord.out.bam  introns_bedfile exons_bedfile output_name_cbc_trimmed_selectedCELLS.txt output_name_cbc_trimmed_star
```
where _introns_bedfile_ and _exons_bedfile_ are bedfiles containing information about the position of introns and exons in the reference genome. 
