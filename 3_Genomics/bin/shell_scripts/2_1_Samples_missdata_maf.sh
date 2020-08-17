#! /bin/sh

#SBATCH -p keri
#SBATCH --mem 80000
#SBATCH -n 2

vcftools --vcf ../../data/TMVB_5SNPradlocus.vcf --keep ../../metadata/88_ind.txt --max-missing 0.9 --maf 0.05 --recode --out ../../data/88ind_maxmiss0.9_maf0.05
