#!/bin/sh

#SBATCH -w nodo3
#SBATCH -n 4

# Veronica Reyes
# Paper:
# Make statistics


for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187$
do samtools flagstat ../../data/BAM/${f}_sw10L50.bam > ../../metadata/Statistics_map/${f}_sw10L50.txt
done
