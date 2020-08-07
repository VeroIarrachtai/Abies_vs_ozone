#!/bin/sh

#SBATCH -w nodo3
#SBATCH -n 4

# Veronica Reyes 
# Paper:
# Make fastQC analyses 


# Do fastqc
for f in ../../data/RAW/*.fastq.gz;
do fastqc --outdir  ../../metadata/fastqc_before_trimm/ $f
done




