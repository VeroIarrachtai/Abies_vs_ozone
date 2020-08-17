#!/bin/sh

#SBATCH -w nodo4
#SBATCH -n 1

# Veronica Reyes
# Paper:
# Make index to alignment


# view the bam file as a regular text file and name it "align":

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do samtools view ../../data/BAM/${f}_sw10L50.bam > ../../data/TXT/SEQ_GENES/${f}_sw10L50.seqgenes.txt
done


# Extract the contig on which each read was aligned on. I think the contig_id is in the third column (please verify with your own bam file before processing, and readjust # -f 3  if necessary). Create a list of contig_id named align1:

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do cut -f 3 ../../data/TXT/SEQ_GENES/${f}_sw10L50.seqgenes.txt > ../../data/TXT/ALL_GENES/${f}_sw10L50.allgenes.txt
done

# Sort the contig list and count how many time each contig_id appears in the list:

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do sort ../../data/TXT/ALL_GENES/${f}_sw10L50.allgenes.txt | uniq -c > ../../data/TXT/COUNT_GENES/${f}_sw10L50.countgenes.txt
done

# Delete "space characters" at the beggining of each line and change the "space" between the number of contig occurrences and the contig_id for a "tab"

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do awk '{ sub(/^[ \t]+/, ""); print }' ../../data/TXT/COUNT_GENES/${f}_sw10L50.countgenes.txt | sed 's/ /\t/' > ../../data/TXT/GENES_ORDER/${f}_sw10L50.genesorder.txt;
done

# There we go, you should obtain a list looking like that:

# 59 AB_000002_T.1
# 23 AB_000003_T.1
# 6 AB_000006_T.1
