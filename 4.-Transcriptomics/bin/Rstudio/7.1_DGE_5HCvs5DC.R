# Verónica Reyes, febrero 2020

# Load libraries
library(limma)
library(edgeR)
library(DESeq2)
library(ggbiplot)
library (ggplot2)

# Load data. Count table 
alldata <-read.delim("../../metadata/all_genes/allreadsgenes.txt")
alldata <- as.data.frame(alldata)

# Convert dataframe to data matrix
x<-alldata
rownames(x)<-alldata[,1] # Add rownames
x<-x[ ,2:ncol(x)] # Remove double col with names
alldata<-as.matrix(x)

########################################### Damaged vs Tolerant 170 ppb###########################################
##################################################################################################################
##################################################################################################################
# Select subset data(descart data)
DCvsHC<- subset(alldata, select = -c(DS_1, DS_2, DS_4,
                                     HS_1, HS_2, HS_5,
                                     HC17, DC47))


############################################
# Add characteristics 
############################################
tratamiento <- c("DC","DC","DC","DC","DC",
                 "HC","HC","HC","HC","HC")
label <- c("DC_1", "DC_2","DC_3","DC_4","DC_5",
          "HC_1","HC_2","HC_3","HC_4","HC_5")
samples <-c("DC1", "DC2","DC3","DC4","DC5",
           "HC1","HC2","HC3","HC4","HC5")
targets <- data.frame(tratamiento,label,samples)
rownames(targets) <- label

targets

### Filtering genes 
table(rowSums(DCvsHC)==0)
suma <- rowSums(DCvsHC)
filtconteos <- DCvsHC[suma != 0,] 
dim(filtconteos)
head(filtconteos)

##################################################################################################################
##################################################################################################################
# EdgeR

## Clase DGEList
d <- DGEList(counts = filtconteos[,1:10], group = targets$tratamiento) ## Normalization
colnames(d) <- targets$label

## Normalization
d <- calcNormFactors(d)
plotMDS(d, main="plotMDS DCvsHC")

## Dispersors stimation
d <- estimateCommonDisp(d,verbose=TRUE)
d <- estimateTagwiseDisp(d)
plotBCV(d, main="plotBCV DCvsHC")

## Test
et <- exactTest(d,pair=c("HC","DC"))
top<- topTags(et, n= Inf)
hist(top$table$FDR, breaks = 100, main = "Hist FDR DCvsHC")
abline(v=0.05, col="red",lwd=3)

##################################################################################################################
##################################################################################################################
# DESeq2

### Class DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=filtconteos, colData= targets, design=~tratamiento)

### Test
dds <- DESeq(dds)
head(dds)
res <- results(dds)
head (res)
### Normalizacion de ambos datos

#edgeR
d$samples$norm.factors #edgeR

#DESeq2
sizeFactors(dds) #DESeq2

####################################
### Compareted dispersion values
####################################
#######
# edgeR
#######

# First calculated comun disspersion
d$common.dispersion
# Second gen to gen dispersion of comun dispersion
head(d$tagwise.dispersion)
# Choose the best of both stimation 

#######
# DESeq2
#######

# First calculated stimation gen to gen  
head(mcols(dds)$dispGeneEs)

# Then through an adjustment with the average counts estimate the dispersion
head(mcols(dds)$dispersion)

############################################################################################################
#Compare the tests, that is, the p-values and other results that each packet has calculated for each gene
############################################################################################################
########
# edgeR
########
head (topTags(et, n= Inf))

########
# DESeq2
########
res[rownames(topTags(et, n= Inf)),]

################################
#Plot Log fold change
################################
########
# edgeR
########p.value=0.1
de <- decideTestsDGE(et, adjust.method = "fdr" )
head(de)
detags <- rownames(d)[as.logical(de)]
plotSmear(et, de.tags=detags, main="plotSmear de edgeR") > abline(h=0, col="red", lwd=3)

########
# DESeq2
plotMA(res, main="MA-plot DESeq2", ylim=c(-5,5))
########

##################################################################
### Sort the genes according to the attached p-value they have obtained
##################################################################
########
# edgeR
########
topSig <- top[top$table$FDR < 0.05, ]
dim(topSig)
head(topSig)
genesDEedgeR <- rownames(topSig)
head(genesDEedgeR)
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
head(topSig_export)
write.table(topSig_export, "../../metadata/DGE/EdgeR_HvsD170ppb_FDR_0.05.txt", sep="\t", row.names=T)

########
# DESeq2
########
# Sort by p-valores
resOrdered <- res[order(res$padj),]
# Only DEG
xx <-res[order(res$padj,na.last=NA),] 
head(xx)
resSig2 <- xx[xx$padj < 0.05, ]
dim(resSig2)
head(resSig2)
genesDEDESeq2 <- rownames(resSig2)
genesDEDESeq2
head(resSig2)
resSig2_export<-resSig2
resSig2_export$ID<-resSig2_export
head(resSig2_export)
write.table(resSig2_export, "../../metadata/DGE/DESeq2_HvsD170ppb_FDR_0.05.txt", sep="\t", row.names=T)

##################################################################
### Sort the genes according to the attached p-value they have obtained
##################################################################
########
# EdgeR
########

topSig <- top[top$table$FDR < 5, ]
dim(topSig)
genesDEedgeR <- rownames(topSig)
head(genesDEedgeR)
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
head(topSig_export)
write.table(topSig_export, "../../metadata/DGE/EdgeR_HvsD170ppb_FDR_5.txt", sep="\t", row.names=T)

########
# DESeq2
########
xx <-res[order(res$padj,na.last=NA),] 
resSig2 <- xx[xx$padj < 5, ]
dim(resSig2)
head(resSig2)
genesDEDESeq2 <- rownames(resSig2)
head(resSig2)
resSig2_export<-resSig2
resSig2_export$ID<-resSig2_export
head(resSig2_export)
write.table(resSig2_export, "../../metadata/DGE/DESeq2_HvsD170ppb_FDR_5.txt", sep="\t", row.names=T)


###################################################
### How many common DE genes exist edgeR vs DESeq2
###################################################

genesDEcomunes <- intersect(genesDEedgeR,genesDEDESeq2) 
head(genesDEcomunes)
str(genesDEcomunes)


