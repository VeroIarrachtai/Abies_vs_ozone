# README Paper Abies vs ozone

This repository contains scripts, data, metadata and results to perform transcriptomic, genetic and metabolic analysis to **"Abies vs ozone's project"**.

Analyzes to answer each particular aim can be found in separate directories with data and pictures:

* Evaluate differential expression of healthy and damaged trees in two ozone's periods (TRANSCRIPTOMICS).

* Identify tolerance sacred fir's origins (GENOMICS).

* Quantify secundary metabolites' relative abundance in healthy and damaged trees during two ozone concenrations' periods (METABOLOMICS).

## Principal directories:

There's repository structure:

```
+----- Abie_vs_ozone/
|	+--1.-Sampling/
|	+--2.-Metabolomics/
|	+--3.-Genomics/
|	+--4.-Transcriptomics/
| +--README.md
```


**README.md**: There's a markdown file about my project. This file includes repository's disposition.

**1.-Sampling**: There's a directory with coordinates and samples' disposition in omics' analisys .

**2.-Metabolomics**: There's metabolites analysis generated with gas chromatograph spectrum mass (GC-SM). Data from html files were loaded into tables. Subsequently relative abundance was calculated. Finally, values between samples were compared from a barplot, ANOVA and PCA analysis.

**3.-Genomics**: There's genomic analysis from GBS sequencing. ipyRAD was used to assemble *de novo*, VCFTools and plink were used to make more specific filters. Relationship was calculated without multiple SNPs in same loci. Mantel test, PCA and admixture were performed and the Heterocity was calculated.

**4.-Transcriptomics**: There's transcriptomic analysis from samples sequenced with RNAseq. Samples were cut with Timmomatic and mapped to a reference transcriptome with BWA. Rstudio allowed to evaluate differential expression between samples with edgeR and DESeq2. Subsequently, volcanoplot was performed to show overexpressed and underexpressed genes.


# Principal analysis workflow:

## GENOMICS

### Identify origins of sacred fir with tolerance to O3.

![](5.-wonderful_images/Genomic_methods.png)

Check more info about this pipeline in the [README_genomics](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/tree/master/1.-GENOMICS/README_genomics.md).

To see a short summary about the analysis of the final data you can go to [Analysis_genomics](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/blob/master/4.-INFO_PROJECT/GENOMICS_ligth_analysis.md)

# METABOLOMICS

### Quantify the relative abundance of secondary metabolites in healthy y damaged trees in two periods of concentration of O3

![](5.-wonderful_images/Metabolomic_methods.png)

Check more info about this pipeline in the [README_metabolomics](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/tree/master/2.-METABOLOMICS/README_metabolomics.md).

To see a short summary about the analysis of the final data you can go to [Analysis_metabolomics](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/blob/master/4.-INFO_PROJECT/METABOLOMICS_ligth_analysis.md)

# TRANSCRIPTOMICS (ALERT: I keep developing this section)

### Evaluate the differential expression of healthy and damaged trees in two periods of [O3].



![](5.-wonderful_images/Transcriptomic_methods.png)

Check more info about this pipeline in the [README_transcriptomics](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/blob/master/3.-TRANSCRIPTOMICS/README_TRANSCRIPTOMICS.md).

To see a short summary about the analysis of the final data you can go to [Analysis_transcriptomics](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/blob/master/4.-INFO_PROJECT/TRANSCRIPTOMICS_ligth_analysis.md)

### Contact
Verónica Reyes Galindo
```
veronica.rg.pb@gmail.com
``
