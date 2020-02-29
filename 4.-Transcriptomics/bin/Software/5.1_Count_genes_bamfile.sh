#!/bin/sh

#SBATCH -p keri
#SBATCH -n 5
#SBATCH --mem=10000


# Verónica Reyes
# 21 Junio 2019
# Contar reads (script from Seb)

# Samples:
# DC01_15_sw10L50_28 DC02_15_sw10L50_28 DC03_15_sw10L50_28 DC04_15_sw10L50_28 DC05_15_sw10L50_28
# DS01_15_sw10L50_28 DS02_15_sw10L50_28 DS04_15_sw10L50_28
# DC04_17_sw10L50_28 SC01_17_sw10L50_28
# SC01_15_sw10L50_28 SC02_15_sw10L50_28 SC03_15_sw10L50_28 SC04_15_sw10L50_28 SC05_15_sw10L50_28
# SS01_15_sw10L50_28 SS02_15_sw10L50_28 SS05_15_sw10L50_28

# Change directory

cd ../../data/BAM


# view the bam file as a regular text file and name it "align":

for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do samtools view ../data/BAM/$i.bam > $i.seqgenes_28.txt; done


# Extract the contig on which each read was aligned on. I think the contig_id is in the third column (please verify with your own bam file before processing, and readjust # -f 3  if necessary). Create a list of contig_id named align1:

for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do cut -f 3 $i.seqgenes_28.txt > $i.allgenes_28.txt; done

# Sort the contig list and count how many time each contig_id appears in the list:

for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do sort $i.allgenes_28.txt | uniq -c > $i.countgenes_28.txt; done


# Delete "space characters" at the beggining of each line and change the "space" between the number of contig occurrences and the contig_id for a "tab"

for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do awk '{ sub(/^[ \t]+/, ""); print }' $i.genesorder_28.txt | sed 's/ /\t/' >; done

# There we go, you should obtain a list looking like that:

# 59 AB_000002_T.1
# 23 AB_000003_T.1
# 6 AB_000006_T.1
# 14 AB_000007_T.1

# Eliminar archivos .bam de Programas

rm *.bam

# Mover todos los archivos creados a su directorio correspondiente

for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/seq_genes/$i.seqgenes_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done
for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/all_genes/$i.allgenes_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done
for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/count_genes/$i.countgenes_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done
for i in DC01_15_sw10L50_28 DC04_17_sw10L50_28 DS04_15_sw10L50_28 SC03_15_sw10L50_28 SS02_15_sw10L50_28 DC02_15_sw10L50_28 DC05_15_sw10L50_28 SC01_15_sw10L50_28 SC04_15_sw10L50_28 SS05_15_sw10L50_28 DC03_15_sw10L50_28 DS01_15_sw10L50_28 SC01_17_sw10L50_28 SC05_15_sw10L50_28 DC04_15_sw10L50_28 DS02_15_sw10L50_28 SC02_15_sw10L50_28 SS01_15_sw10L50_28; do mv ../../metadata/genes_order/$i.genesorder_28.txt ../../TRANSCRIPTOMICS_MAP/Count/Map_REFTransSNP; done
