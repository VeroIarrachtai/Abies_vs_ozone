# README Paper Abies vs ozone

This repository contains scripts, data, metadata, and results to perform transcriptomic, genetic, and metabolic analysis to **"Abies vs ozone's project"**.

In these directories you will find analysis to answer these project's particular aims:

* To evaluate the healthy and damaged trees' differential expression in two ozone's periods (TRANSCRIPTOMICS).

* To identify tolerance Sacred fir's origins (GENOMICS).

* To quantify secundary metabolites' relative abundance in healthy and damaged trees during two ozone concenrations' periods (METABOLOMICS).

## Principal directories:

There's repository structure:

```
+----- Abies_vs_ozone/
|	+--1.-Sampling/
|	+--2.-Metabolomics/
|	+--3.-Genomics/
|	+--4.-Transcriptomics/
|	+--5.-INFO_PROJECT/
|	+--README.md
```


**README.md**: There's a markdown file about my project. This file includes repository's disposition.

**1.-Sampling**: There's a directory with coordinates and samples' disposition in omics' analisys .

**2.-Metabolomics**: There're metabolites analysis generated with gas chromatograph spectrum mass (GC-SM). Data from html files were loaded into tables. Subsequently relative abundance was calculated. Finally, values between samples were compared from a barplot, ANOVA and PCA analysis.

**3.-Genomics**: There're genomic analysis from GBS sequencing. ipyRAD was used to assemble *de novo*, VCFTools and plink were used to make more specific filters. Relationship was calculated without multiple SNPs in same loci. Mantel test, PCA and admixture were performed and the Heterocity was calculated.

**4.-Transcriptomics**: There're transcriptomic analysis from samples sequenced with RNAseq. Samples were cut with Timmomatic and mapped to a reference transcriptome with BWA. Rstudio allowed to evaluate differential expression between samples with edgeR and DESeq2. Subsequently, volcanoplot was performed to show overexpressed and underexpressed genes.

**5.-INFO_PROJECT**: There're exhibitions, summaries and final analysis.

# Principal analysis workflow:

## GENOMICS

### Identify tolerance sacred fir's origins (GENOMICS).

![](3.-Genomics/metadata/Genomic_methods.png)


Check more info about this pipeline in the [README_genomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/3.-Genomics/README_genomics.md).

To see a short summary about the analysis of the final data you can to go here: [Analysis_genomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/5.-INFO_PROJECT/GENOMICS_ligth_analysis.md)

# METABOLOMICS

### Quantify secundary metabolites' relative abundance in healthy and damaged trees during two ozone concenrations' periods (METABOLOMICS).

![](2.-Metabolomics/metadata/Metabolomic_methods.png)

Check more info about this pipeline in the [README_metabolomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/2.-Metabolomics/README_metabolomics.md).

To see a short summary about the analysis of the final data you can to go here: [Analysis_metabolomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/5.-INFO_PROJECT/METABOLOMICS_ligth_analysis.md)

### Evaluate differential expression of healthy and damaged trees in two ozone's periods (TRANSCRIPTOMICS).

![](4.-Transcriptomics/metadata/Transcriptomic_methods.png)

Check more info about this pipeline in the [README_transcriptomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/4.-Transcriptomics/README_TRANSCRIPTOMICS.md).

To see a short summary about the analysis of the final data you can to go here: [Analysis_transcriptomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/5.-INFO_PROJECT/TRANSCRIPTOMICS_ligth_analysis.md)

### Contact

```
Verónica Reyes Galindo
```
