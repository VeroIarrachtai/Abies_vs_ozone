#!/bin/sh

#SBATCH -w nodo3
#SBATCH -n 4

# Veronica Reyes
# Paper:
# Make index to alignment


# To view bam files as a regular text file we are using "view". Later we saved this view like "txt" files. Final files named "*.seqgenes.txt":

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do samtools view ../../data/BAM/${f}_sw10L50.bam > ../../metadata/seq_genes/${f}_sw10L50.seqgenes.txt
done


# To extract the contig on which each read was aligned on. I think the contig_id is in the third column (please verify with your own bam file before processing, and readjust # -f 3  if necessary). Create a list of contig_id named align1:

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do cut -f 3 ../../metadata/seq_genes/${f}_sw10L50.seqgenes.txt > ../../metadata/all_genes/${f}_sw10L50.allgenes.txt
done

# Sort the contig list and count how many time each contig_id appears in the list:

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do sort ../../metadata/all_genes/${f}_sw10L50.allgenes.txt | uniq -c > ../../metadata/count_genes/${f}_sw10L50.countgenes.txt
done

cd ../../metadata/count_genes

# Delete "space characters" at the beggining of each line and change the "space" between the number of contig occurrences and the contig_id for a "tab"

for f in DPVR1_S179 DPVR2_S180 DPVR3_S181 DPVR4_S182 DPVR5_S183 DPVR6_S184 DPVR7_S185 DPVR8_S186 DPVR9_S187 DPVR10_S188 DPVR11_S189 DPVR12_S190 DPVR13_S191 DPVR14_S192 DPVR15_S193 DPVR16_S194 DPVR17_S195 DPVR18_S196;
do awk '{ sub(/^[ \t]+/, ""); print }' ${f}_sw10L50.genesorder.txt | sed 's/ /\t/' >;
done

mv *genes_order.txt ../../metadata/genes_order/

cd ../../bin/Software

# There we go, you should obtain a list looking like that:

# 59 AB_000002_T.1
# 23 AB_000003_T.1
# 6 AB_000006_T.1
