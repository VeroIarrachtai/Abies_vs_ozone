#!/bin/sh

#SBATCH -w nodo3
#SBATCH -n 4

# Veronica Reyes 
# Paper: Abies vs ozone
# Make index to alignment

bwa index -p ../../metadata/Index/index_Areligiosa -a is ../../metadata/Reference_Transcriptome/GCAT_AB-RNA-1.0.16.fa 
