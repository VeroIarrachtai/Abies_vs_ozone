#!/bin/sh

#SBATCH -w nodo3
#SBATCH -n 4

# Veronica Reyes
# Paper:
# Make statistics


#Create list with names over and down

# Extract data BEST ORF

for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_C_TvsD_IDs_over.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_C_TvsD_over.txt
for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_C_TvsD_IDs_down.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_C_TvsD_down.txt

for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_S_TvsD_IDs_over.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_S_TvsD_over.txt
for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_S_TvsD_IDs_down.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_S_TvsD_down.txt

for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_D_170Cvs87SS_IDs_over.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_D_170Cvs87SS_over.txt
for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_D_170Cvs87SS_IDs_down.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_D_170Cvs87SS_down.txt

for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_T_170Cvs87SS_IDs_over.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_T_170Cvs87SS_over.txt
for i in $(cat ../../data/Over_Down/SPECIFIC/L_D2_ER_T_170Cvs87SS_IDs_down.txt) ;do grep -A1 "$i" ../../../../../../Descargas/GCAT_AB-RNA-1.0.16.Pfam.txt; done > ../../data/Over_Down/SPECIFIC/P-fam_D2_ER_T_170Cvs87SS_down.txt
