# README Paper Abies vs ozone

This repository contains the scripts, data, metadata, and results to perform the transcriptomic, the genetic, and the metabolic analysis to **"Abies vs Ozone's Project"**.

In these directories you will find analyses to answer the project's particular aims:

* To quantify the secondary metabolite relative abundance in the tolerant and the damaged trees during two ozone concentration periods: [Metabolomic analysis](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/2_Metabolomics)

* To identify the tolerance sacred fir's origins:  [Genomic analysis](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/3_Genomics)

* To evaluate the tolerant and the damaged trees differential expression in the two ozone concentration periods: [Transcriptomic analysis](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/4_Transcriptomics)


## GENERAL directory structure:

```
+----- Abies_vs_ozone/
|	+--1_Sampling/
|	+--2_Metabolomics/
|	+--3_Genomics/
|	+--4_Transcriptomics/
|	+--5_INFO_PROJECT/
|	+--README.md
```

**README.md**: a markdown file about this project. This file includes the repository's disposition.

**1_Sampling**: a directory with the coordinates and the samples' disposition in the omics analisys.

**2_Metabolomics**: metabolite-generated analyses with a gas chromatograph spectrum mass (GC-SM). Data from the html files were loaded into the tables. Subsequently, the relative abundance was calculated. Finally, all the samples' values were compared using a barplot, ANOVA, and PCA analysis.

**3_Genomics**: genomic analyses from the GBS sequencing. **ipyRAD** was used to assemble *de novo*, **VCFTools** and **plink** were used to make more specific filters. The relatedness was calculated using a single SNP per locus. The mantel test, PCA, and admixture were performed in order to identify samples' local origin.

**4_Transcriptomics**: transcriptomic analysis from RNAseq data. Samples were cut with **Trimmomatic** and they were mapped to a reference transcriptome with **BWA**. **R** allowed to evaluate differential expression between the samples with the packages **edgeR** and **DESeq2**. Subsequently, a volcanoplot was performed to show overexpressed and underexpressed genes.

**5_INFO_PROJECT**: slide-shows, summaries, and final analyses.

# Workfolow of main analyses:

## [METABOLOMICS](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/2_Metabolomics)

### To quantify the secondary metabolite relative abundance in the healthy and the damaged trees during the two ozone concenrations periods.

![](2_Metabolomics/metadata/Metabolomic_methods.png)

Check more information about this pipeline in [README_metabolomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/2_Metabolomics/README.md).

To see a short summary about the final metabolomics analysis click [here](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/5_INFO_PROJECT/METABOLOMICS_ligth_analysis.md).

## [GENOMICS](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/3_Genomics)

### To identify the geographic origin of tolerant invididuals.

![](3_Genomics/metadata/Genomic_methods.png)

Check more information about this pipeline in [README_genomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/3_Genomics/README.md).

To see a short summary about the final genomic analysis click [here](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/5_INFO_PROJECT/GENOMICS_ligth_analysis.md).

## [TRANSCRIPTOMICS](https://github.com/VeroIarrachtai/Abies_vs_ozone/tree/master/4_Transcriptomics)

### To evaluate the healthy and the damaged trees differential expression in the two periods with contrasting ozone concentrations.

![](4_Transcriptomics/metadata/Transcriptomic_methods.png)

Check more information about this pipeline in [README_transcriptomics](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/4_Transcriptomics/README.md).

To see a short summary about the final transcriptomic analysis click [here](https://github.com/VeroIarrachtai/Abies_vs_ozone/blob/master/5_INFO_PROJECT/TRANSCRIPTOMICS_ligth_analysis.md).

### Contact

```
Verónica Reyes Galindo
veronica.rg.pb@gmail.com
```
