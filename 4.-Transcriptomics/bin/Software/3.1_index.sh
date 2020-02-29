#!/bin/bash

#SBATCH -p keri
#SBATCH -n 1
#SBATCH --mem=150000

bwa index -p ../../metadata/INDEX/index_Areligiosa -a is ../../metadata/Reference_Transcriptome/GCAT_AB-RNA-1.0.16.fa
