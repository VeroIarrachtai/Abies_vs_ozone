#! /bin/sh

#SBATCH -p keri
#SBATCH --mem 80000
#SBATCH -n 2

./plink --bfile ../../data/snp_withoutDupLoci_without_duplicates89s_maxmiss0.9_maf0.05 --make-rel square --make-bed --out ../../data/relsnp_withoutDupLoci_without_duplicates88s_maxmiss0.9_maf0.05
./plink --bfile ../../data/relsnp_withoutDupLoci_without_duplicates89s_maxmiss0.9_maf0.05 --recode --out ../../data/relsnp_withoutDupLoci_without_duplicates88s_maxmiss0.9_maf0.05
./plink --file ../../data/relsnp_withoutDupLoci_without_duplicates89s_maxmiss0.9_maf0.05 --recode vcf --out ../../data/relsnp_withoutDupLoci_without_duplicates88s_maxmiss0.9_maf0.05
