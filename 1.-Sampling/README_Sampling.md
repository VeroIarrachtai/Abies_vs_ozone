# README SAMPLING

Before starting the analysis here are the programs that need to be installed:

Before to start the analysis you need to install these programs:

## Prerequisites

VCFtools 0.1.15
PLINK v1.9
R 3.4.2
Admixture 1.3

### R packages:
SNPStats
ggplot2
SNPRelate

# Repository structure:

## Principal directories:
```
+---- Abies_religiosa_vs_ozone/
|	+--1.-Sampling/
|          +--bin/
|          +--data/
|          +--metadata/
|          +--outputs/
|          +--README_Sampling.md
```

**README.md**: is a intro about my project. This include the structure of this repository.

**1.-Sampling**

**bin**:

**data**: It is a genomic analysis from samples product of a GBS sequencing, ipyRAD was used to assemble de novo, VCFTools and plink to make more specific filters. The relationship was calculated, multiple SNPs were discarded in the same loci, a Mantel test, PCA and admixture were performed and the Heterocity was calculated.

**metadata**: This is a transcriptomic analysis from samples sequenced with RNAseq. Samples were cut(Timmomatic) and mapped (BWA) to a reference transcriptome. Sequence counting was carried out through command lines in Rstudio that subsequently allowed the evaluation of differential expression between samples. From the counting table, a volcanoplot was performed to exemplify the overexpressed and underexpressed genes.

**outputs**


### Contact
```
Verónica Reyes Galindo
veronica.rg.pb@gmail.com
```


Data S2: specific number of readings and the percentage of mapping per individual and condition
