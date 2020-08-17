# Reyes Galindo Verónica
# 12 Febrero 2019

### This script draw a plot with eigenvects values using SNPRelate


## Load packages
library(gdsfmt)
library(SNPRelate)
library(ggplot2)

#Load data file (".vcf")
vcf.fn <- "../../data/without_Dup_loci/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05.vcf"

## Reformat to ".gds" file
snpgdsVCF2GDS(vcf.fn, "../../outputs/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05_pca.gds", method="biallelic.only", verbose = TRUE)

## Get summary ".gds" file
snpgdsSummary("../../outputs/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05_pca.gds")

## Open ".gds" file
genofile <- snpgdsOpen("../../outputs/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05_pca.gds")

## Run PCA
pca<-snpgdsPCA(genofile,remove.monosnp=TRUE, num.thread=2)

## Get variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

## Make a data.frame with eigenvects values
tab <- data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1],    # the first eigenvector,
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

## Load metadata file
fullmat<-read.csv("../../metadata/PLACA_FINAL_89_samples.csv", header = TRUE, sep = ",")

## Match metadata info with eigenvectors values
tab$Poblacion<-fullmat$Localidad[match(tab$sample.id, fullmat$key_comun)]

## Choose nice colors and ozone condition
pobcol<- c( "#e94b57",
            "#000000",
            "#76c74f")
tab$Condicion<- c(rep("otro",7),rep("damaged",4),
                  rep("tolerant",5), rep("otro",72))

## Draw PCA with ozone condition
ggplot(tab, aes(x=EV1, y=EV2))+
  geom_point(aes(color=Condicion, shape=Poblacion), size =5) +
  scale_color_manual(values = pobcol) +
  theme(legend.title = element_text(size=15))+
  theme(legend.text = element_text(size = 15))+
  xlab(paste0("Eigenvector 1 explaining ", round(pc.percent, 2)[1], "%")) +
  ylab(paste0("Eigenvector 2 explaining ", round(pc.percent, 2)[2], "%"))+
  theme(axis.title.y = element_text(size = rel(2), angle = 90))+
  theme(axis.title.x = element_text(size = rel(2), angle = 360))+
  theme(axis.text.x = element_text(hjust = .5, size=13, color="black"))+
  theme(axis.text.y = element_text(hjust = .5, size=13, color="black"))+
  geom_point(alpha = 1/20)+
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,11,13,15,16,17,18,35,38,43,64))

ggsave("../../outputs/6.1_PCA.png")
