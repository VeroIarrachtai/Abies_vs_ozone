#! /bin/sh

#SBATCH -p keri
#SBATCH --mem 80000
#SBATCH -n 2

vcftools --vcf ../../data/snp_withoutDupLoci_89s_maxmiss0.9_maf0.05.vcf --keep ../../metadata/samples_het_relat.txt --het --out ../../data/samples_he_snp_withoutDupLoci_10ind_maxmiss0.9_maf0.05.het
