# Reyes Galindo Verónica
# 12 Febrero 2019

### This script calculated IBD

## Load packages
library(geosphere)
library(gdsfmt) 
library(SNPRelate)
library(ggplot2)

## Load packages
require(MASS) 
library(vegan)
library(permute)
library(lattice)
library(psych)

#Load metadata
Coordenadas<-read.delim("../../metadata/Ar_IBD2.txt")

# Convert dataframe to matrix
matcoor<- as.data.frame(Coordenadas)
row.names(matcoor)
x<-matcoor
rownames(x)<-matcoor[,1]
x<-x[ ,2:ncol(x)]
x
rownames(x)
matcoor<-as.matrix(x)

## Calculate geographic distances

#Sierra Manantlan vs Todas
matcoor<-matcoor[-1,]
SM<-distGeo(matcoor, c(-103.9500, 19.4500), a=6378137, f=1/298.257223563)
# Aj vs Todas
matcoor<-matcoor[-1,]
Aj<-distGeo(matcoor, c(-99.2340, 19.2230), a=6378137, f=1/298.257223563)
# Santa Rosa Xochiac vs Todas 
matcoor<-matcoor[-1,]
SRX<-distGeo(matcoor, c(-99.3010, 19.2850), a=6378137, f=1/298.257223563)
# El Chico vs Todas 
matcoor<-matcoor[-1,]
EC<-distGeo(matcoor, c(-98.7000, 20.1500), a=6378137, f=1/298.257223563)
# Nevado de Colima vs Todas 
matcoor<-matcoor[-1,]
NC<-distGeo(matcoor, c(-103.5990, 19.5840), a=6378137, f=1/298.257223563)
# Nevado de Toluca NT vs Todas 
matcoor<-matcoor[-1,]
NT<-distGeo(matcoor, c(-99.8100, 19.1870), a=6378137, f=1/298.257223563)
# Ixtapalucan vs Todas 
matcoor<-matcoor[-1,]
Ix<-distGeo(matcoor, c(-98.6100, 19.2590), a=6378137, f=1/298.257223563)
# Cerro Blanco vs Todas 
matcoor<-matcoor[-1,]
CB<-distGeo(c(-100.2280, 19.5650),matcoor, a=6378137, f=1/298.257223563)
# Puerta Garnica vs Todas 
matcoor<-matcoor[-1,]
PG<-distGeo(matcoor, c(-100.8220, 19.6700), a=6378137, f=1/298.257223563)
# Volcan Tancitaro   vs Todas 
matcoor<-matcoor[-1,]
VT<-distGeo(matcoor, c(-102.3170, 19.3830), a=6378137, f=1/298.257223563)
# Michoacan vs Todas 
matcoor<-matcoor[-1,]
MiAl<-distGeo(matcoor, c(-100.6040, 19.8000), a=6378137, f=1/298.257223563)
# San Andres vs Todas 
matcoor<-matcoor[-1,]
SA<-distGeo(matcoor, c(-100.6040, 19.8000), a=6378137, f=1/298.257223563)
# Nevado de Toluca RG vs Todas 
matcoor<-matcoor[-1,]
NTRG<-distGeo(matcoor, c(-99.9600, 19.2600), a=6378137, f=1/298.257223563)
# Nevado de Toluca SB vs Todas 
matcoor<-matcoor[-1,]
NTSB<-distGeo(matcoor, c(-99.9200, 19.2300), a=6378137, f=1/298.257223563)
# Volcan Atlitzin  vs Todas 
matcoor<-matcoor[-1,]
VA<-distGeo(matcoor, c(-97.3500, 18.9670), a=6378137, f=1/298.257223563)
# Cerro Zamorano vs Todas 
matcoor<-matcoor[-1,]
CZ<-distGeo(matcoor, c(-100.1820, 20.9290), a=6378137, f=1/298.257223563)
# Tlaxco  vs Todas 
matcoor<-matcoor[-1,]
Tl<-distGeo(matcoor, c(-98.0830, 19.6830), a=6378137, f=1/298.257223563)
# Malinche  vs Todas 
matcoor<-matcoor[-1,]
Ma<-distGeo(matcoor, c(-98.0433, 19.2372), a=6378137, f=1/298.257223563)
#Cofre de Perote vs Todas
#matcoor<-matcoor[-1,]
#CP<-distGeo(matcoor, c(-97.1500, 19.5170), a=6378137, f=1/298.257223563)
allvectordist<- c(SM,Aj, SRX, EC, NC, NT,Ix, CB, PG, VT, MiAl,
                  SA, NTRG, NTSB, VA, CZ, Tl, Ma)


## Make ".gds" file
vcf.fn <- "../../data/without_Dup_loci/snp_withoutDupLoci_89ind_maxmiss0.9_maf0.05.vcf"
snpgdsVCF2GDS(vcf.fn, "../../outputs/snp_withoutDupLoci_89ind_maxmiss0.9_maf0.05.recode.gds", method="biallelic.only", verbose = TRUE)
snpgdsSummary("../../outputs/snp_withoutDupLoci_89ind_maxmiss0.9_maf0.05.recode.gds")
genofile<-snpgdsOpen("../../outputs/snp_withoutDupLoci_89ind_maxmiss0.9_maf0.05.recode.gds")

# Load population data 
pop_code <- read.delim("../../metadata/FST_VCFTools_Ar89.txt", header=FALSE)


# Get sample id
sample.id<- read.gdsn(index.gdsn(genofile, "sample.id"))

## Estimate FST
##############
##"SM", "Aj"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "Aj") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsAj<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsAj$Fst
##############
##"SM", "SRX"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "SRX") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsSRX<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsSRX$Fst
##############
##"SM", "EC"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "EC") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsEC<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsEC$Fst
##############
##"SM", "NC"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "NC") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsNC<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsNC$Fst
##############
##"SM", "NT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "NT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsNT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsNT$Fst
##############
##"SM", "Ix"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "Ix") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsIx<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsIx$Fst
##############
##"SM", "CB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "CB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsCB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsCB$Fst
##############
##"SM", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsPG$Fst
##############
##"Aj", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsVT$Fst
##############
##"SM", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsMiAl$Fst
##############
##"SM", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsSA$Fst
##############
##"SM", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsNTRG$Fst
##############
##"SM", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsNTSB$Fst
##############
##"SM", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsVA$Fst
##############
##"SM", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsCZ$Fst
##############
##"SM", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsTl$Fst
##############
##"SM", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsMa$Fst
##############
##"SM", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SM", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SMvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SMvsCP$Fst



##############
##"Aj", "SRX"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "SRX") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsSRX<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsSRX$Fst
##############
##"Aj", "EC"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "EC") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsEC<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsEC$Fst
##############
##"Aj", "NC"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "NC") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsNC<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsNC$Fst
##############
##"Aj", "NT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "NT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsNT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsNT$Fst
##############
##"Aj", "Ix"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "Ix") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsIx<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsIx$Fst
##############
##"Aj", "CB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "CB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsCB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsCB$Fst
##############
##"Aj", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsPG$Fst
##############
##"Aj", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsVT$Fst
##############
##"Aj", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsMiAl$Fst
##############
##"Aj", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsSA$Fst
##############
##"Aj", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsNTRG$Fst
##############
##"Aj", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsNTSB$Fst
##############
##"Aj", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsVA$Fst
##############
##"Aj", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsCZ$Fst
##############
##"Aj", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsTl$Fst
##############
##"Aj", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsMa$Fst
##############
##"Aj", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Aj", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
AjvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
AjvsCP$Fst

#SANTA ROSA XOCHIAC

##############
##"SRX", "EC"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "EC") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsEC<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsEC$Fst
##############
##"SRX", "NC"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "NC") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsNC<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsNC$Fst
##############
##"SRX", "NT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "NT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsNT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsNT$Fst
##############
##"SRX", "Ix"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "Ix") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsIx<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsIx$Fst
##############
##"SRX", "CB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "CB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsCB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsCB$Fst
##############
##"SRX", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsPG$Fst
##############
##"SRX", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsVT$Fst
##############
##"SRX", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsMiAl$Fst
##############
##"SRX", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsSA$Fst
##############
##"SRX", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsNTRG$Fst
##############
##"SRX", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsNTSB$Fst
##############
##"SRX", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsVA$Fst
##############
##"SRX", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsCZ$Fst
##############
##"SRX", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsTl$Fst
##############
##"SRX", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsMa$Fst
##############
##"SRX", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SRX", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SRXvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SRXvsCP$Fst

#EL CHICO

##############
##"EC", "NC"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "NC") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsNC<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsNC$Fst
##############
##"EC", "NT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "NT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsNT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsNT$Fst
##############
##"EC", "Ix"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "Ix") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsIx<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsIx$Fst
##############
##"EC", "CB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "CB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsCB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsCB$Fst
##############
##"EC", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsPG$Fst
##############
##"EC", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsVT$Fst
##############
##"EC", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsMiAl$Fst
##############
##"EC", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsSA$Fst
##############
##"EC", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsNTRG$Fst
##############
##"EC", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsNTSB$Fst
##############
##"EC", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsVA$Fst
##############
##"EC", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsCZ$Fst
##############
##"EC", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsTl$Fst
##############
##"EC", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsMa$Fst
##############
##"EC", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("EC", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
ECvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
ECvsCP$Fst


# Nevado de Colima

##############
##"NC", "NT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "NT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsNT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsNT$Fst
##############
##"NC", "Ix"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "Ix") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsIx<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsIx$Fst
##############
##"NC", "CB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "CB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsCB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsCB$Fst
##############
##"NC", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsPG$Fst
##############
##"NC", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsVT$Fst
##############
##"NC", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsMiAl$Fst
##############
##"NC", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsSA$Fst
##############
##"NC", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsNTRG$Fst
##############
##"NC", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsNTSB$Fst
##############
##"NC", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsVA$Fst
##############
##"NC", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsCZ$Fst
##############
##"NC", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsTl$Fst
##############
##"NC", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsMa$Fst
##############
##"NC", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NC", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NCvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NCvsCP$Fst


# NEVADO DE TOLUCA


##############
##"NT", "Ix"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "Ix") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsIx<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsIx$Fst
##############
##"NT", "CB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "CB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsCB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsCB$Fst
##############
##"NT", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsPG$Fst
##############
##"NT", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsVT$Fst
##############
##"NT", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsMiAl$Fst
##############
##"NT", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsSA$Fst
##############
##"NT", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsNTRG$Fst
##############
##"NT", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsNTSB$Fst
##############
##"NT", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsVA$Fst
##############
##"NT", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsCZ$Fst
##############
##"NT", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsTl$Fst
##############
##"NT", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsMa$Fst
##############
##"NT", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NT", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTvsCP$Fst


#IXTAPALUCAN

##############
##"Ix", "CB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "CB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsCB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsCB$Fst
##############
##"Ix", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsPG$Fst
##############
##"Ix", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsVT$Fst
##############
##"Ix", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsMiAl$Fst
##############
##"Ix", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsSA$Fst
##############
##"Ix", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsNTRG$Fst
##############
##"Ix", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsNTSB$Fst
##############
##"Ix", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsVA$Fst
##############
##"Ix", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsCZ$Fst
##############
##"Ix", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsTl$Fst
##############
##"Ix", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsMa$Fst
##############
##"Ix", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ix", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
IxvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
IxvsCP$Fst

#CERRO BLANCO

##############
##"CB", "PG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "PG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsPG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsPG$Fst
##############
##"CB", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsVT$Fst
##############
##"CB", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsMiAl$Fst
##############
##"CB", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsSA$Fst
##############
##"CB", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsNTRG$Fst
##############
##"CB", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsNTSB$Fst
##############
##"CB", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsVA$Fst
##############
##"CB", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsCZ$Fst
##############
##"CB", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsTl$Fst
##############
##"CB", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsMa$Fst
##############
##"CB", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CB", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CBvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CBvsCP$Fst

# PUERTA GARNICA

##############
##"PG", "VT"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "VT") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsVT<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsVT$Fst
##############
##"PG", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsMiAl$Fst
##############
##"PG", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsSA$Fst
##############
##"PG", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsNTRG$Fst
##############
##"PG", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsNTSB$Fst
##############
##"PG", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsVA$Fst
##############
##"PG", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsCZ$Fst
##############
##"PG", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsTl$Fst
##############
##"PG", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsMa$Fst
##############
##"PG", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("PG", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
PGvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
PGvsCP$Fst

# VOLCAN TANCITARO


##############
##"VT", "MiAl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "MiAl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsMiAl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsMiAl$Fst
##############
##"VT", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsSA$Fst
##############
##"VT", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsNTRG$Fst
##############
##"VT", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsNTSB$Fst
##############
##"VT", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsVA$Fst
##############
##"VT", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsCZ$Fst
##############
##"VT", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsTl$Fst
##############
##"VT", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsMa$Fst
##############
##"VT", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VT", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VTvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VTvsCP$Fst

# MICHOACAN

##############
##"MiAl", "SA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "SA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsSA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsSA$Fst
##############
##"MiAl", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsNTRG$Fst
##############
##"MiAl", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsNTSB$Fst
##############
##"MiAl", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsVA$Fst
##############
##"MiAl", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsCZ$Fst
##############
##"MiAl", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsTl$Fst
##############
##"MiAl", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsMa$Fst
##############
##"MiAl", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("MiAl", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MiAlvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MiAlvsCP$Fst

# SAN ANDRES

##############
##"SA", "NTRG"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SA", "NTRG") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SAvsNTRG<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SAvsNTRG$Fst
##############
##"SA", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SA", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SAvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SAvsNTSB$Fst
##############
##"SA", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SA", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SAvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SAvsVA$Fst
##############
##"SA", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SA", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SAvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SAvsCZ$Fst
##############
##"SA", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SA", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SAvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SAvsTl$Fst
##############
##"SA", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SA", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SAvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SAvsMa$Fst
##############
##"SA", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("SA", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
SAvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
SAvsCP$Fst

# NEVADO DE TOLUCA RG

##############
##"NTRG", "NTSB"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTRG", "NTSB") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTRGvsNTSB<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTRGvsNTSB$Fst
##############
##"NTRG", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTRG", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTRGvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTRGvsVA$Fst
##############
##"NTRG", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTRG", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTRGvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTRGvsCZ$Fst
##############
##"NTRG", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTRG", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTRGvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTRGvsTl$Fst
##############
##"NTRG", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTRG", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTRGvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTRGvsMa$Fst
##############
##"NTRG", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTRG", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTRGvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTRGvsCP$Fst

# NEVADO DE TOLUCA SB

##############
##"NTSB", "VA"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTSB", "VA") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTSBvsVA<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTSBvsVA$Fst
##############
##"NTSB", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTSB", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTSBvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTSBvsCZ$Fst
##############
##"NTSB", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTSB", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTSBvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTSBvsTl$Fst
##############
##"NTSB", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTSB", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTSBvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTSBvsMa$Fst
##############
##"NTSB", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("NTSB", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
NTSBvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
NTSBvsCP$Fst

# VOLCAN ATLITZIN

##############
##"VA", "CZ"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VA", "CZ") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VAvsCZ<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VAvsCZ$Fst
##############
##"VA", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VA", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VAvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VAvsTl$Fst
##############
##"VA", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VA", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VAvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VAvsMa$Fst
##############
##"VA", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("VA", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
VAvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
VAvsCP$Fst

#CERRO ZAMORANO

##############
##"CZ", "Tl"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CZ", "Tl") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CZvsTl<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CZvsTl$Fst
##############
##"CZ", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CZ", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CZvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CZvsMa$Fst
##############
##"CZ", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("CZ", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
CZvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
CZvsCP$Fst

#TLAXCO

##############
##"Tl", "Ma"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Tl", "Ma") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
TlvsMa<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
TlvsMa$Fst
##############
##"Tl", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Tl", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
TlvsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
TlvsCP$Fst

# MALINCHE

##############
##"Ma", "CP"##
##############
# Flag desired populations
flag <-  pop_code$V2 %in% c("Ma", "CP") 
samp.sel <- sample.id[flag]
pop.sel <- pop_code$V2[flag]
pop.sel<-droplevels(pop.sel)
# Estimate FST
MavsCP<-snpgdsFst(genofile, sample.id=samp.sel, population= as.factor(pop.sel), autosome.only = FALSE, remove.monosnp=TRUE, method="W&C84")
MavsCP$Fst

## Make matrix with FST data
allvectorfst<- c(SMvsSRX$Fst,SMvsSRX$Fst,SMvsEC$Fst,SMvsNC$Fst,SMvsNT$Fst,SMvsIx$Fst,SMvsCB$Fst,SMvsPG$Fst,SMvsVT$Fst,SMvsMiAl$Fst,SMvsSA$Fst,SMvsNTRG$Fst,SMvsNTSB$Fst,SMvsVA$Fst,SMvsCZ$Fst,SMvsTl$Fst,SMvsMa$Fst,SMvsCP$Fst,
                 AjvsSRX$Fst,AjvsEC$Fst,AjvsNC$Fst,AjvsNT$Fst,AjvsIx$Fst,AjvsCB$Fst,AjvsPG$Fst,AjvsVT$Fst,AjvsMiAl$Fst,AjvsSA$Fst,AjvsNTRG$Fst,AjvsNTSB$Fst,AjvsVA$Fst,AjvsCZ$Fst,AjvsTl$Fst,AjvsMa$Fst,AjvsCP$Fst,
                 SRXvsEC$Fst,SRXvsNC$Fst,SRXvsNT$Fst,SRXvsIx$Fst,SRXvsCB$Fst,SRXvsPG$Fst,SRXvsVT$Fst,SRXvsMiAl$Fst,SRXvsSA$Fst,SRXvsNTRG$Fst,SRXvsNTSB$Fst,SRXvsVA$Fst,SRXvsCZ$Fst,SRXvsTl$Fst,SRXvsMa$Fst,SRXvsCP$Fst,
                 ECvsNC$Fst,ECvsNT$Fst,ECvsIx$Fst,ECvsCB$Fst,ECvsPG$Fst,ECvsVT$Fst,ECvsMiAl$Fst,ECvsSA$Fst,ECvsNTRG$Fst,ECvsNTSB$Fst,ECvsVA$Fst,ECvsCZ$Fst,ECvsTl$Fst,ECvsMa$Fst,ECvsCP$Fst,
                 NCvsNT$Fst,NCvsIx$Fst,NCvsCB$Fst,NCvsPG$Fst,NCvsVT$Fst,NCvsMiAl$Fst,NCvsSA$Fst,NCvsNTRG$Fst,NCvsNTSB$Fst,NCvsVA$Fst,NCvsCZ$Fst,NCvsTl$Fst,NCvsMa$Fst,NCvsCP$Fst,
                 NTvsIx$Fst,NTvsCB$Fst,NTvsPG$Fst,NTvsVT$Fst,NTvsMiAl$Fst,NTvsSA$Fst,NTvsNTRG$Fst,NTvsNTSB$Fst,NTvsVA$Fst,NTvsCZ$Fst,NTvsTl$Fst,NTvsMa$Fst,NTvsCP$Fst,
                 IxvsCB$Fst,IxvsPG$Fst,IxvsVT$Fst,IxvsMiAl$Fst,IxvsSA$Fst,IxvsNTRG$Fst,IxvsNTSB$Fst,IxvsVA$Fst,IxvsCZ$Fst,IxvsTl$Fst,IxvsMa$Fst,IxvsCP$Fst,
                 CBvsPG$Fst,CBvsVT$Fst,CBvsMiAl$Fst,CBvsSA$Fst,CBvsNTRG$Fst,CBvsNTSB$Fst,CBvsVA$Fst,CBvsCZ$Fst,CBvsTl$Fst,CBvsMa$Fst,CBvsCP$Fst,
                 PGvsVT$Fst,PGvsMiAl$Fst,PGvsSA$Fst,PGvsNTRG$Fst,PGvsNTSB$Fst,PGvsVA$Fst,PGvsCZ$Fst,PGvsTl$Fst,PGvsMa$Fst,PGvsCP$Fst,
                 VTvsMiAl$Fst,VTvsSA$Fst,VTvsNTRG$Fst,VTvsNTSB$Fst,VTvsVA$Fst,VTvsCZ$Fst,VTvsTl$Fst,VTvsMa$Fst,VTvsCP$Fst,
                 MiAlvsSA$Fst,MiAlvsNTRG$Fst,MiAlvsNTSB$Fst,MiAlvsVA$Fst,MiAlvsCZ$Fst,MiAlvsTl$Fst,MiAlvsMa$Fst,MiAlvsCP$Fst,
                 SAvsNTRG$Fst,SAvsNTSB$Fst,SAvsVA$Fst,SAvsCZ$Fst,SAvsTl$Fst,SAvsMa$Fst,SAvsCP$Fst,
                 NTRGvsNTSB$Fst,NTRGvsVA$Fst,NTRGvsCZ$Fst,NTRGvsTl$Fst,NTRGvsMa$Fst,NTRGvsCP$Fst,
                 NTSBvsVA$Fst,NTSBvsCZ$Fst,NTSBvsTl$Fst,NTSBvsMa$Fst,NTSBvsCP$Fst,
                 VAvsCZ$Fst,VAvsTl$Fst,VAvsMa$Fst,VAvsCP$Fst,
                 CZvsTl$Fst,CZvsMa$Fst,CZvsCP$Fst,
                 TlvsMa$Fst,TlvsCP$Fst,
                 MavsCP$Fst
)

IBD<-read.delim("../../metadata/Ar_IBD_comparations.txt")
IBD$FST<-allvectorfst
IBD$Dist<-allvectordist

# Add Distance in Km.

IBD$Dist_km <- IBD$Dist/1000

## Draw IBD
ggplot(data = IBD, aes(x = IBD$Dist_km, y = IBD$FST)) + 
  geom_point(colour = "#e60000") +
  xlab("Distancia geográfica (Km)")+
  ylab("Distancia Genética")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_smooth(method=lm)+
  theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.5), angle = 360))+
  theme(axis.text.x = element_text(hjust = .5, size=10, color="black"))+
  theme(axis.text.y = element_text(hjust = .5, size=10, color="black"))

## Make  mantel test
pairs(x = IBD, lower.panel = NULL)

## Load packages

pairs.panels(x = IBD, ellipses = FALSE, lm = TRUE, method = "pearson")

## Make a matrix to mantel test
SierraManantlan<-                                                                                                                                                                                                                               c(0,SMvsAj$Fst,SMvsSRX$Fst,SMvsEC$Fst,SMvsNC$Fst,SMvsNT$Fst,SMvsIx$Fst,SMvsCB$Fst,SMvsPG$Fst,SMvsVT$Fst,SMvsMiAl$Fst,SMvsSA$Fst,SMvsNTRG$Fst,SMvsNTSB$Fst,SMvsVA$Fst,SMvsCZ$Fst,SMvsTl$Fst,SMvsMa$Fst,SMvsCP$Fst)
Ajusco<-                                                                                                                                                                                                                               c(SMvsAj$Fst,0,AjvsSRX$Fst,AjvsEC$Fst,AjvsNC$Fst,AjvsNT$Fst,AjvsIx$Fst,AjvsCB$Fst,AjvsPG$Fst,AjvsVT$Fst,AjvsMiAl$Fst,AjvsSA$Fst,AjvsNTRG$Fst,AjvsNTSB$Fst,AjvsVA$Fst,AjvsCZ$Fst,AjvsTl$Fst,AjvsMa$Fst,AjvsCP$Fst)
SantaRosaXochiac<-                                                                                                                                                                                                           c(SMvsAj$Fst,AjvsSRX$Fst,0,SRXvsEC$Fst,SRXvsNC$Fst,SRXvsNT$Fst,SRXvsIx$Fst,SRXvsCB$Fst,SRXvsPG$Fst,SRXvsVT$Fst,SRXvsMiAl$Fst,SRXvsSA$Fst,SRXvsNTRG$Fst,SRXvsNTSB$Fst,SRXvsVA$Fst,SRXvsCZ$Fst,SRXvsTl$Fst,SRXvsMa$Fst,SRXvsCP$Fst)
ElChico<-                                                                                                                                                                                                           c(SMvsAj$Fst,AjvsEC$Fst,SRXvsEC$Fst,0,ECvsNC$Fst,ECvsNT$Fst,ECvsIx$Fst,ECvsCB$Fst,ECvsPG$Fst,ECvsVT$Fst,ECvsMiAl$Fst,ECvsSA$Fst,ECvsNTRG$Fst,ECvsNTSB$Fst,ECvsVA$Fst,ECvsCZ$Fst,ECvsTl$Fst,ECvsMa$Fst,ECvsCP$Fst) 
NevadoColima<-                                                                                                                                                                                             c(SMvsAj$Fst,AjvsNC$Fst,SRXvsNC$Fst,ECvsNC$Fst,0,NCvsNT$Fst,NCvsIx$Fst,NCvsCB$Fst,NCvsPG$Fst,NCvsVT$Fst,NCvsMiAl$Fst,NCvsSA$Fst,NCvsNTRG$Fst,NCvsNTSB$Fst,NCvsVA$Fst,NCvsCZ$Fst,NCvsTl$Fst,NCvsMa$Fst,NCvsCP$Fst)
NevadoToluca<-                                                                                                                                                                                    c(SMvsAj$Fst,AjvsNT$Fst,SRXvsNT$Fst,ECvsNT$Fst,NCvsNT$Fst,0,NTvsIx$Fst,NTvsCB$Fst,NTvsPG$Fst,NTvsVT$Fst,NTvsMiAl$Fst,NTvsSA$Fst,NTvsNTRG$Fst,NTvsNTSB$Fst,NTvsVA$Fst,NTvsCZ$Fst,NTvsTl$Fst,NTvsMa$Fst,NTvsCP$Fst)
Ixtapaluca<-                                                                                                                                                                            c(SMvsAj$Fst,AjvsIx$Fst,SRXvsIx$Fst,ECvsIx$Fst,NCvsIx$Fst,NTvsIx$Fst,0,IxvsCB$Fst,IxvsPG$Fst,IxvsVT$Fst,IxvsMiAl$Fst,IxvsSA$Fst,IxvsNTRG$Fst,IxvsNTSB$Fst,IxvsVA$Fst,IxvsCZ$Fst,IxvsTl$Fst,IxvsMa$Fst,IxvsCP$Fst)
CerroBlanco<-                                                                                                                                                                 c(SMvsAj$Fst,AjvsCB$Fst,SRXvsCB$Fst,ECvsCB$Fst,NCvsCB$Fst,NTvsCB$Fst,IxvsCB$Fst,0,CBvsPG$Fst,CBvsVT$Fst,CBvsMiAl$Fst,CBvsSA$Fst,CBvsNTRG$Fst,CBvsNTSB$Fst,CBvsVA$Fst,CBvsCZ$Fst,CBvsTl$Fst,CBvsMa$Fst,CBvsCP$Fst)
PuertaGarnica<-                                                                                                                                                      c(SMvsAj$Fst,AjvsPG$Fst,SRXvsPG$Fst,ECvsPG$Fst,NCvsPG$Fst,NTvsPG$Fst,IxvsPG$Fst,CBvsPG$Fst,0,PGvsVT$Fst,PGvsMiAl$Fst,PGvsSA$Fst,PGvsNTRG$Fst,PGvsNTSB$Fst,PGvsVA$Fst,PGvsCZ$Fst,PGvsTl$Fst,PGvsMa$Fst,PGvsCP$Fst)
VolcanTancitaro<-                                                                                                                                           c(SMvsAj$Fst,AjvsVT$Fst,SRXvsVT$Fst,ECvsVT$Fst,NCvsVT$Fst,NTvsVT$Fst,IxvsVT$Fst,CBvsVT$Fst,PGvsVT$Fst,0,VTvsMiAl$Fst,VTvsSA$Fst,VTvsNTRG$Fst,VTvsNTSB$Fst,VTvsVA$Fst,VTvsCZ$Fst,VTvsTl$Fst,VTvsMa$Fst,VTvsCP$Fst)
Michoacan<-                                                                                                                      c(SMvsAj$Fst,AjvsMiAl$Fst,SRXvsMiAl$Fst,ECvsMiAl$Fst,NCvsMiAl$Fst,NTvsMiAl$Fst,IxvsMiAl$Fst,CBvsMiAl$Fst,PGvsMiAl$Fst,VTvsMiAl$Fst,0,MiAlvsSA$Fst,MiAlvsNTRG$Fst,MiAlvsNTSB$Fst,MiAlvsVA$Fst,MiAlvsCZ$Fst,MiAlvsTl$Fst,MiAlvsMa$Fst,MiAlvsCP$Fst)
SanAndres<-                                                                                                                            c(SMvsAj$Fst,AjvsSA$Fst,SRXvsSA$Fst,ECvsSA$Fst,NCvsSA$Fst,NTvsSA$Fst,IxvsSA$Fst,CBvsSA$Fst,PGvsSA$Fst,VTvsSA$Fst,MiAlvsSA$Fst,0,SAvsNTRG$Fst,SAvsNTSB$Fst,SAvsVA$Fst,SAvsCZ$Fst,SAvsTl$Fst,SAvsMa$Fst,SAvsCP$Fst)
NevadoTolucaRG<-                                                                                       c(SMvsAj$Fst,AjvsNTRG$Fst,SRXvsNTRG$Fst,ECvsNTRG$Fst,NCvsNTRG$Fst,NTvsNTRG$Fst,IxvsNTRG$Fst,CBvsNTRG$Fst,PGvsNTRG$Fst,VTvsNTRG$Fst,MiAlvsNTRG$Fst,SAvsNTRG$Fst,0,NTRGvsNTSB$Fst,NTRGvsVA$Fst,NTRGvsCZ$Fst,NTRGvsTl$Fst,NTRGvsMa$Fst,NTRGvsCP$Fst)
NevadoTolucaSB<-                                                                          c(SMvsAj$Fst,AjvsNTSB$Fst,SRXvsNTSB$Fst,ECvsNTSB$Fst,NCvsNTSB$Fst,NTvsNTSB$Fst,IxvsNTSB$Fst,CBvsNTSB$Fst,PGvsNTSB$Fst,VTvsNTSB$Fst,MiAlvsNTSB$Fst,SAvsNTSB$Fst,NTRGvsNTSB$Fst,0,NTSBvsVA$Fst,NTSBvsCZ$Fst,NTSBvsTl$Fst,NTSBvsMa$Fst,NTSBvsCP$Fst)
VolcanAtlitzin<-                                                                                       c(SMvsAj$Fst,AjvsVA$Fst,SRXvsVA$Fst,ECvsVA$Fst,NCvsVA$Fst,NTvsVA$Fst,IxvsVA$Fst,CBvsVA$Fst,PGvsVA$Fst,VTvsVA$Fst,MiAlvsVA$Fst,SAvsVA$Fst,NTRGvsVA$Fst,NTSBvsVA$Fst,0,VAvsCZ$Fst,VAvsTl$Fst,VAvsMa$Fst,VAvsCP$Fst)
CerroZamorano<-                                                                               c(SMvsAj$Fst,AjvsCZ$Fst,SRXvsCZ$Fst,ECvsCZ$Fst,NCvsCZ$Fst,NTvsCZ$Fst,IxvsCZ$Fst,CBvsCZ$Fst,PGvsCZ$Fst,VTvsCZ$Fst,MiAlvsCZ$Fst,SAvsCZ$Fst,NTRGvsCZ$Fst,NTSBvsCZ$Fst,VAvsCZ$Fst,0,CZvsTl$Fst,CZvsMa$Fst,CZvsCP$Fst)
Tlaxco<-                                                                             c(SMvsAj$Fst,AjvsTl$Fst,SRXvsTl$Fst,ECvsTl$Fst,NCvsTl$Fst,NTvsTl$Fst,IxvsTl$Fst,CBvsTl$Fst,PGvsTl$Fst,VTvsTl$Fst,MiAlvsTl$Fst,SAvsTl$Fst,NTRGvsTl$Fst,NTSBvsTl$Fst,VAvsTl$Fst,CZvsTl$Fst,0,TlvsMa$Fst,TlvsCP$Fst)
Malinche<-                                                                  c(SMvsAj$Fst,AjvsMa$Fst,SRXvsMa$Fst,ECvsMa$Fst,NCvsMa$Fst,NTvsMa$Fst,IxvsMa$Fst,CBvsMa$Fst,PGvsMa$Fst,VTvsMa$Fst,MiAlvsMa$Fst,SAvsMa$Fst,NTRGvsMa$Fst,NTSBvsMa$Fst,VAvsMa$Fst,CZvsMa$Fst,TlvsMa$Fst,0,MavsCP$Fst)
CofrePerote<-                                                      c(SMvsAj$Fst,AjvsCP$Fst,SRXvsCP$Fst,ECvsCP$Fst,NCvsCP$Fst,NTvsCP$Fst,IxvsCP$Fst,CBvsCP$Fst,PGvsCP$Fst,VTvsCP$Fst,MiAlvsCP$Fst,SAvsCP$Fst,NTRGvsCP$Fst,NTSBvsCP$Fst,VAvsCP$Fst,CZvsCP$Fst,TlvsCP$Fst,MavsCP$Fst,0)

mat <- cbind(SierraManantlan,Ajusco, SantaRosaXochiac, ElChico,NevadoColima,NevadoToluca,
             Ixtapaluca,CerroBlanco,PuertaGarnica,VolcanTancitaro,Michoacan,
             SanAndres,NevadoTolucaRG,NevadoTolucaSB,VolcanAtlitzin,
             CerroZamorano,Tlaxco, Malinche,CofrePerote
)

mat2<- mat/(1-mat)

Fst.dists<- as.matrix(mat)
mountain.dists <- dist(cbind(x$Lon, x$Lat))
mountain.dists <- as.matrix(mountain.dists)[1:19, 1:19]

Fst.dists<- as.dist(mat)
mountain.dists <-as.dist(mountain.dists)

mantel(mountain.dists, Fst.dists, method="pearson", permutations=999)

########################################################################################################################################################################

# Linearize as suggested by Rousset (1997) for IBD using FST/(1 − FST)
IBD$FST_liner<- IBD$FST/(1-IBD$FST)

# Add Distance in Km.

IBD$Dist_km <- IBD$Dist/1000


## Draw IBD using FST/(1 − FST)
ggplot(data = IBD, aes(x = IBD$Dist_km, y = IBD$FST_liner)) + 
  geom_point(colour = "black") +
  xlab("Geographic distance (Km)")+
  ylab("Genetic distance")+
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_smooth(method=lm)+
  theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.5), angle = 360))+
  theme(axis.text.x = element_text(hjust = .5, size=10, color="black"))+
  theme(axis.text.y = element_text(hjust = .5, size=10, color="black"))

ggsave("../../outputs/5.1_Mantel_test_Linear.png")

## Make  mantel test
pairs(x = IBD, lower.panel = NULL)

## Load packages

pairs.panels(x = IBD, ellipses = FALSE, lm = TRUE, method = "pearson")

## Make a matrix to mantel test
mat2<- mat/(1-mat)

Fst.dists<- as.matrix(mat2)
mountain.dists <- dist(cbind(x$Lon, x$Lat))
mountain.dists <- as.matrix(mountain.dists)[1:19, 1:19]

Fst.dists<- as.dist(mat2)
mountain.dists <-as.dist(mountain.dists)

mantel(mountain.dists, Fst.dists, method="pearson", permutations=999)





