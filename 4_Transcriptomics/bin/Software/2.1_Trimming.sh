#!/bin/sh

#SBATCH -w nodo3
#SBATCH -n 4

# Veronica Reyes 
# Paper:
# Make trimming and fastQC analyses 


# Do Trimmer  with Trimmomatic-0.36

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do java -jar ../../../../Programas/Trimmomatic/Trimmomatic_bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33 ../../data/RAW/${f}_L007_R1_001.fastq.gz ../../data/RAW/${f}_L007_R2_001.fastq.gz ../../data/TRIMMING/Trimmer_${f}_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimmer_${f}_L007_R1_001_unpaired.fq.gz ../../data/TRIMMING/Trimmer_${f}_L007_R2_001_paired.fq.gz ../../data/TRIMMING/Trimmer_${f}_L007_R2_001_unpaired.fq.gz ILLUMINACLIP:../../../../Programas/Trimmomatic/Trimmomatic_bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:10:28 MINLEN:50 HEADCROP:13
done

# Do fastqc after timming

for f in ../../data/TRIMMING/*.fq.gz;
do fastqc --outdir  ../../metadata/fastqc_after_trimm/ $f
done
