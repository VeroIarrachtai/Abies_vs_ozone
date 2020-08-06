#! /bin/sh

#SBATCH -p keri
#SBATCH --mem 80000
#SBATCH -n 2

vcftools --vcf ../../data/89ind_maxmiss0.9_maf0.05.recode.vcf --freq --out ../../data/freq_89ind_maxmiss0.9_maf0.05
