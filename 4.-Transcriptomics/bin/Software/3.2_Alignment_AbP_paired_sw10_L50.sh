#!/bin/sh

#SBATCH -w nodo3
#SBATCH -n 6

# Veronica Reyes
# Paper:
# Make alignment

for a in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do bwa mem ../../metadata/Index/index_Areligiosa ../../data/TRIMMING/Trimmer_${a}_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimmer_${a}_L007_R2_001_paired.fq.gz > ../../data/SAM/${a}_sw10L50.sam
done
