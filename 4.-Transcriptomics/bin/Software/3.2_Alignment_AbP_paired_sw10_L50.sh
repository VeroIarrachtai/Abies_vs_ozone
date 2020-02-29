#!/bin/bash

#SBATCH -p keri
#SBATCH -n 5
#SBATCH --mem=500000

bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR1_S179_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR1_S179_L007_R2_001_paired.fq.gz > ../../data/SAM/SC01_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR2_S180_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR2_S180_L007_R2_001_paired.fq.gz > ../../data/SAM/SC02_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR3_S181_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR3_S181_L007_R2_001_paired.fq.gz > ../../data/SAM/SC03_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR4_S182_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR4_S182_L007_R2_001_paired.fq.gz > ../../data/SAM/SC04_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR5_S183_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR5_S183_L007_R2_001_paired.fq.gz > ../../data/SAM/SC05_15_sw10L50_28.sam

bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR6_S184_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR6_S184_L007_R2_001_paired.fq.gz > ../../data/SAM/DC01_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR7_S185_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR7_S185_L007_R2_001_paired.fq.gz > ../../data/SAM/DC02_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR8_S186_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR8_S186_L007_R2_001_paired.fq.gz > ../../data/SAM/DC03_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR9_S187_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR9_S187_L007_R2_001_paired.fq.gz > ../../data/SAM/DC04_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR10_S188_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR10_S188_L007_R2_001_paired.fq.gz > ../../data/SAM/DC05_15_sw10L50_28.sam

bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR11_S189_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR11_S189_L007_R2_001_paired.fq.gz > ../../data/SAM/SS01_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR12_S190_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR12_S190_L007_R2_001_paired.fq.gz > ../../data/SAM/SS02_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR13_S191_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR13_S191_L007_R2_001_paired.fq.gz > ../../data/SAM/SS05_15_sw10L50_28.sam

bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR14_S192_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR14_S192_L007_R2_001_paired.fq.gz > ../../data/SAM/DS01_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR15_S193_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR15_S193_L007_R2_001_paired.fq.gz > ../../data/SAM/DS02_15_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR16_S194_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR16_S194_L007_R2_001_paired.fq.gz > ../../data/SAM/DS04_15_sw10L50_28.sam

bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR17_S195_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR17_S195_L007_R2_001_paired.fq.gz > ../../data/SAM/SC01_17_sw10L50_28.sam
bwa mem ../../metadata/index_GCAT_AB-RNA-1.0.16/index_Areligiosa ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR18_S196_L007_R1_001_paired.fq.gz ../../data/TRIMMING/Trimm18s_sw10-28_ml50_28/Trimmer_DPVR18_S196_L007_R2_001_paired.fq.gz > ../../data/SAM/DC04_17_sw10L50_28.sam
