# Verónica Reyes, febrero 2020
# Damaged vs Tolerant 170 ppb

# Load libraries
library(limma)
library(edgeR)
library(DESeq2)
library(ggbiplot)
library (ggplot2)

# Load count table 
alldata <-read.delim("../../metadata/all_genes/allreadsgenes.txt")

# Convert dataframe to data matrix
x<-alldata
rownames(x)<-alldata[,1] # Add rownames
x<-x[ ,2:ncol(x)] # Remove double col with names
alldata<-as.matrix(x)

# Select subset data(descart data)
DCvsTC<- subset(alldata, select = -c(DS_1, DS_2, DS_4,
                                     TS_1, TS_2, TS_5,
                                     TC17, DC47))

# Add characteristics 

treatment <- c("DC","DC","DC","DC","DC",
                 "TC","TC","TC","TC","TC")
label <- c("DC_1", "DC_2","DC_3","DC_4","DC_5",
          "TC_1","TC_2","TC_3","TC_4","TC_5")
samples <-c("DC1", "DC2","DC3","DC4","DC5",
           "TC1","TC2","TC3","TC4","TC5")
targets <- data.frame(treatment,label,samples)
rownames(targets) <- label
targets

# Filtering genes 
table(rowSums(DCvsTC)==0)
suma <- rowSums(DCvsTC)
filtconteos <- DCvsTC[suma != 0,] 
dim(filtconteos)
head(filtconteos)

##################################################################################################################
# EdgeR
## Clase DGEList
d <- DGEList(counts = filtconteos[,1:10], group = targets$treatment) ## Normalization
colnames(d) <- targets$label

## Normalization
d <- calcNormFactors(d)
plotMDS(d, main="plotMDS DCvsTC")

## Dispersors stimation
d <- estimateCommonDisp(d,verbose=TRUE)
d <- estimateTagwiseDisp(d)
plotBCV(d, main="plotBCV DCvsTC")

## Test
et <- exactTest(d,pair=c("TC","DC"))
top<- topTags(et, n= Inf)
hist(top$table$FDR, breaks = 100, main = "Hist FDR DCvsTC")
abline(v=0.05, col="red",lwd=3)

##################################################################################################################
# DESeq2
### Class DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=filtconteos, colData= targets, design=~treatment)

### Test
dds <- DESeq(dds)
head(dds)
res <- results(dds)
head (res)

### Normalizacion de ambos datos
# edgeR
d$samples$norm.factors #edgeR
# DESeq2
sizeFactors(dds) #DESeq2

### Compareted dispersion values
# edgeR
# First calculated comun disspersion
d$common.dispersion
# Second gen to gen dispersion of comun dispersion
head(d$tagwise.dispersion)
# Choose the best of both stimation 

# DESeq2
# First calculated stimation gen to gen  
head(mcols(dds)$dispGeneEs)
# Then through an adjustment with the average counts estimate the dispersion
head(mcols(dds)$dispersion)

#Compare the tests, that is, the p-values and other results that each packet has calculated for each gene
# edgeR
head (topTags(et, n= Inf))
# DESeq2
res[rownames(topTags(et, n= Inf)),]


#Plot Log fold change
# edgeR
########p.value=0.1
de <- decideTestsDGE(et, adjust.method = "fdr" )
head(de)
detags <- rownames(d)[as.logical(de)]
plotSmear(et, de.tags=detags, main="plotSmear de edgeR") > abline(h=0, col="red", lwd=3)

# DESeq2
plotMA(res, main="MA-plot DESeq2", ylim=c(-5,5))

### Sort the genes according to the attached p-value they have obtained
# edgeR
topSig <- top[top$table$FDR < 0.0, ]
dim(topSig)
head(topSig)
genesDEedgeR <- rownames(topSig)
head(genesDEedgeR)
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
head(topSig_export)
write.table(topSig_export, "../../metadata/DGE/EdgeR_HvsD170ppb_FDR_0.01.txt", sep="\t", row.names=T)

# DESeq2
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

### Sort the genes according to the attached p-value they have obtained
# EdgeR
topSig <- top[top$table$FDR < 5, ]
dim(topSig)
genesDEedgeR <- rownames(topSig)
head(genesDEedgeR)
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
head(topSig_export)
write.table(topSig_export, "../../metadata/DGE/EdgeR_HvsD170ppb_FDR_5.txt", sep="\t", row.names=T)

# DESeq2
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

### How many common DE genes exist edgeR vs DESeq2
genesDEcomunes <- intersect(genesDEedgeR,genesDEDESeq2) 
head(genesDEcomunes)
str(genesDEcomunes)

