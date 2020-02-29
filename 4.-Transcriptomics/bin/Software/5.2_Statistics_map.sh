#!/bin/sh

#SBATCH -p keri
#SBATCH -n 5
#SBATCH --mem=10000


# Verónica Reyes
# 21 Junio 2019
# Contar reads (script from Seb)

./samtools flagstat ../../data/BAM/DC01_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DC01_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/DC02_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DC02_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/DC03_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DC03_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/DC04_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DC04_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/DC05_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DC05_15_sw10L50_TR.txt

./samtools flagstat ../../data/BAM/DC04_17_sw10L50_TR.bam > ../../metadata/Statistics_map/DC04_17_sw10L50_TR.txt

./samtools flagstat ../../data/BAM/SC01_17_sw10L50_TR.bam > ../../metadata/Statistics_map/SC01_17_sw10L50_TR.txt

./samtools flagstat ../../data/BAM/SC01_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SC01_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/SC02_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SC02_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/SC03_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SC03_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/SC04_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SC04_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/SC05_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SC05_15_sw10L50_TR.txt

./samtools flagstat ../../data/BAM/DS01_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DS01_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/DS02_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DS02_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/DS04_15_sw10L50_TR.bam > ../../metadata/Statistics_map/DS04_15_sw10L50_TR.txt

./samtools flagstat ../../data/BAM/SS01_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SS01_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/SS02_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SS02_15_sw10L50_TR.txt
./samtools flagstat ../../data/BAM/SS05_15_sw10L50_TR.bam > ../../metadata/Statistics_map/SS05_15_sw10L50_TR.txt
