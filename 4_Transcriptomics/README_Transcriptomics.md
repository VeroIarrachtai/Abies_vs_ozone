# README TRANSCRIPTOMICS

## Pre-requisitos

Before to start the genomic analysis you will need to install:

## SOFTWARE

* [bwa](http://bio-bwa.sourceforge.net)
* [samtools](http://www.htslib.org)
* [Trimmomatic-0.36](http://www.usadellab.org/cms/?page=trimmomatic)
* [R](https://cran.r-project.org)
* [Rstudio (optional)](https://rstudio.com)

## R packages

* **VennDiagram**
* **limma**
* **edgeR**
* **DESeq2**

## TRASNCRIPTOMICS directory structure:

```
+----- Abies_religiosa_vs_ozone/
|	+--3_TRASNCRIPTOMICS/
|		+--bin/
|	   	    +--Rstudio/
|	   	       +--6_1_Countreads_makematrix.R
|	   	       +--7_1_5HCvs5DC.R
|	   	       +--8_1_Volcanoplot.R
|	   	    +--Software/
|	   	       +--1_1_FastQC.sh
|	   	       +--2_1_Trimming.sh
|	   	       +--3_1_index.sh
|	   	       +--3_2_Alignment_AbP_paired_sw10_L50.sh
|	   	       +--4_1_ConvertSamBam_sw10_L50.sh
|	   	       +--5_1_Count_genes_bamfile.sh
|	   	       +--5_2_Statistics_map.sh
|		+--data/
|	   	    +--BAM/
|	   	       +--*.bam
|	   	    +--RAW/
|	   	       +--DPVR1_S179_L007_R1_001_fastqc.zip
|	   	    +--SAM/
|	   	       +--*.sam
|	   	    +--TRIMMING/
|	   	       +--Trimm18s_sw10-28_ml50_28/
|	   	           +--Trimmer_DPVR1_S179_L007_R1_001_paired_fastqc.zip
|	   	    +--TXT/
|	   	       +--*.allgenes.txt
|	   	       +--*.countgenes.txt
|	   	       +--*.genesorder.txt
|	   	       +--*.seqgenes.txt
|	   	    +--OSF.md
|	   	    +--DGE/
|	   	       +--DESeq2_HvsD170ppb_FDR_0.05.txt
|	   	       +--EdgeR_HvsD170ppb_FDR_0.05.txt
|		+--metadata/
|	   	    +--Statistics_map/
|	   	       +--flagstat_DC01_15_sw10L50.txt
|	   	    +--fastqc_before_trimm/
|	   	       +--DPVR1_S179_L007_R1_001_fastqc.html
|	   	    +--fastqc_after_trimm/
|	   	       +--DPVR1_S179_L007_R1_001_fastqc.html
|	   	    +--genes_order/
|	   	       +--DC01_15_sw10L50.genesorder.txt
|	   	    +--Index/
|	   	       +--index_GCAT_AB-RNA-1.0.16/
|	   	       +--index_Areligiosa
|	   	    +--Reference_Trancriptome/
|	   	       +--GCAT_AB-RNA-1.0.16.fa
|	   	    +--RNA_sacred_fir.csv
|	   	    +--laneBarcode.html
|		+--outputs/
|		+--README_transcriptomics.md
|		+--Transcriptomic_methods.png
```

## TRANSCRIPTOMICS content

:file_folder: **`/bin`**
There are scripts and data necessary to do transcriptomic analysis.
Inside of this directory there are two directories. **Software** directory contains scripts that run in [terminal](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/tree/master/4.-Transcriptomics/bin/Software) using command line. **Rstudio** directory contains the scripts that run in [Rstudio](https://github.com/VeroIarrachtai/Abies_religiosa_vs_ozone/tree/master/4.-Transcriptomics/bin/Rstudio).

:file_folder: **`/data`**
There are the producted files of the sequencing and analysis of them (bam, sam, and fastqc.zip files).

:file_folder: **`/metadata`** There are tables and data that complement the omics data. Such as name of samples, count per genes, name of genes, name of sequences, index, reference transcriptome, etc.

:file_folder: **`/outputs`** The figures from Rstudio are stored here.

:page_facing_up: **`/README_transcriptomics`** This is a README that describes the steps to perform the data analysis. It is organized numerically. It is explained that input is necessary and what outputs are obtained from each step.

# 1.0.- Evaluate sequences with fastqc

* **INPUT**:
   * **fitered_file.vcf**(88ind_maxmiss0.9_maf0.05.recode.vcf)

* **OUTPUT**:
   * **fitered_file.freq**(freq_88ind_maxmiss0.9_maf0.05.frq)
   * **fitered_file.bed**(88ind_maxmiss0.9_maf0.05.bed)

## 1.1.-You have to get how often the loci have

SCRIPT in 4_Transcriptomics/Software/[1_1_FastQC.sh](bin/Software/1_1_FastQC.sh)

```


```

**OUT: barplot_images.png**

# 2.0.- Cut sequences with Trimmomatic

* **INPUT**:
   * **fitered_file.vcf**(88ind_maxmiss0.9_maf0.05.recode.vcf)

* **OUTPUT**:
   * **fitered_file.freq**(freq_88ind_maxmiss0.9_maf0.05.frq)
   * **fitered_file.bed**(88ind_maxmiss0.9_maf0.05.bed)

## 2.1.-

SCRIPT in 4_Transcriptomics/Software/[2_1_Trimming.sh](bin/Software/2_1_Trimming.sh)

```
#Do Trimmer  with Trimmomatic-0.36

java -jar ../../Programas/Trimmomatic/Trimmomatic_bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 ../data/RAW/DPVR1_S179_L007_R1_001.fastq.gz ../data/RAW/DPVR1_S179_L007_R2_001.fastq.gz Trimmer_DPVR1_S179_L007_R1_001_paired.fq.gz Trimmer_DPVR1_S179_L007_R1_001_unpaired.fq.gz Trimmer_DPVR1_S179_L007_R2_001_paired.fq.gz Trimmer_DPVR1_S179_L007_R2_001_unpaired.fq.gz ILLUMINACLIP:../../Programas/Trimmomatic/Trimmomatic_bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:10:28 MINLEN:50 HEADCROP:13

```

**OUT:Trimmer_fastqc.png**

# 3.0.- BWA mapping

* **INPUT**:
   * **fitered_file.vcf**(88ind_maxmiss0.9_maf0.05.recode.vcf)

* **OUTPUT**:
   * **fitered_file.freq**(freq_88ind_maxmiss0.9_maf0.05.frq)
   * **fitered_file.bed**(88ind_maxmiss0.9_maf0.05.bed)

## 3.1.- Do Index

SCRIPT in 4_Transcriptomics/Software/[3_1_index.sh](bin/Software/3_1_index.sh)

```
bwa index -p ../../metadata/INDEX/index_Areligiosa -a is ../../metadata/Reference_Transcriptome/GCAT_AB-RNA-1.0.16.fa
```
**OUT: index_Areligiosa.amb, .ann, .bwt, .pac, .sa**

## 3.2.- Mapeo en A. balsamea

SCRIPT in 4_Transcriptomics/Software/[3_2_Alignment_AbP_paired_sw10_L50](bin/Software/3_2_Alignment_AbP_paired_sw10_L50.sh)

```
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR1_S179_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR1_S179_L007_R2_001_paired.fq.gz > ../../data/SAM/SC01_15_sw10L50_28.sam
```
**OUT: files.sam**

# 4.0.- Convert SAM to BAM

* **INPUT**:
   * **fitered_file.vcf**(88ind_maxmiss0.9_maf0.05.recode.vcf)

* **OUTPUT**:
   * **fitered_file.freq**(freq_88ind_maxmiss0.9_maf0.05.frq)
   * **fitered_file.bed**(88ind_maxmiss0.9_maf0.05.bed)

## 4.1.-

SCRIPT in 4_Transcriptomics/Software/[4_1_ConvertSamBam_sw10_L50](bin/Software/4_1_ConvertSamBam_sw10_L50.sh)

```
cd ../../data/sam

for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do samtools view -Sb ../../data/SAM/$i.sam > ../../data/BAM/$i.bam; done

cd ../../bin/Software

```

**OUT: files.bam**

# 5.0.- Count genes in bam file

* **INPUT**:
   * **fitered_file.vcf**(88ind_maxmiss0.9_maf0.05.recode.vcf)

* **OUTPUT**:
   * **fitered_file.freq**(freq_88ind_maxmiss0.9_maf0.05.frq)
   * **fitered_file.bed**(88ind_maxmiss0.9_maf0.05.bed)

## 5.1.-

SCRIPT in 4_Transcriptomics/Software/[5_1_Count_genes_bamfile](bin/Software/5_1_Count_genes_bamfile.sh)

```
cd ../../data/BAM

for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/seq_genes/$i.seqgenes_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done
for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/all_genes/$i.allgenes_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done
for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/count_genes/$i.countgenes_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done
for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/genes_order/$i.genesorder_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done

```

**OUT:**

## 5.2.-Statistics con BWA


SCRIPT in 4_Transcriptomics/Software/[5_2_Statistics_map](bin/Software/5_2_Statistics_map.sh)

```
./samtools flagstat ../../data/BAM/DC01_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DC01_15_sw10L50_TR.txt
```

# 6.0.- Table of transcript counts

* **INPUT**:
   * **genesorder.txt**(DC01_15_sw10L50.genesorder.txt, DC02_15_sw10L50.genesorder.txt, DC03_15_sw10L50.genesorder.txt, etc.)

* **OUTPUT**:
   * **allreadsgenes.txt**(allreadsgenes.txt)

## 6.1.- Table of transcript counts

SCRIPT in 4_Transcriptomics/Rstudio/[6_1_Countreads_makematrix.R](bin/Rstudio/6_1_Countreads_makematrix.R)

**OUT: allreadsgenes.txt**

![count_Table](../5_wonderful_images/Count_table.png)

# 7.0.- DGE analysis

* **INPUT**:
  * **allreadsgenes.txt**(allreadsgenes.txt)

* **OUTPUT**:
  * **images.png**(images.png)

## 7.1.- DGE analysis

SCRIPT in 4_Transcriptomics/Rstudio/[7_1_5HCvs5DC.R](bin/Rstudio/7_1_5HCvs5DC.R)

**OUT: DESeq2.txt, EdgeR.txt**

![](../5_wonderful_images/DESeq2.png)


# 8.0.- Volcano plot

* **INPUT**:
  * **allreadsgenes.txt**(allreadsgenes.txt)

* **OUTPUT**:
  * **images.png**(images.png)

## 8.1.- Volcano plot

SCRIPT in 4_Transcriptomics/Rstudio/[8_1_Volcanoplot.R](bin/Rstudio/8_1_Volcanoplot.R)

**OUT: images.png**

![](outputs/8_1_VPSol_DESeq2.png)
![](outputs/8_1_VPSol_edge.png)

#Contacto

```
Verónica Reyes Galindo
veronica.rg.pb@gmail.com
```
