#! /bin/sh

#SBATCH -p keri
#SBATCH --mem 80000
#SBATCH -n 2

./plink --file ../../data/89ind_maxmiss0.9_maf0.05 --extract ../../metadata/positions_s89_Ar0.9.txt  --make-bed --out ../../data/snp_withoutDupLoci_89s_maxmiss0.9_maf0.05
