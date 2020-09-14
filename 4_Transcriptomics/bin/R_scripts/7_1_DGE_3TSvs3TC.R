# Verónica Reyes, febrero 2020

# Load libraries
library(limma)
library(edgeR)
library(DESeq2)
library(ggbiplot)
library (ggplot2)

# Load data. Count table 
alldata <-read.delim("../../data/allreadsgenes.txt")
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
TSvsTC<- subset(alldata, select = c(TS_1, TS_2, TS_5,
                                     TC_1, TC_2, TC_4))


############################################
# Add characteristics 
############################################
tratamiento <- c("TS","TS","TS",
                 "TC","TC","TC")
label <- c("TS_1", "TS_2", "TS_5",
           "TC_1", "TC_2", "TC_4")
samples <-c("TS1", "TS2", "TS5",
            "TC1", "TC2", "TC4")
targets <- data.frame(tratamiento,label,samples)
rownames(targets) <- label

targets

### Filtering genes 
table(rowSums(TSvsTC)==0)
suma <- rowSums(TSvsTC)
filtconteos <- TSvsTC[suma != 0,] 
dim(filtconteos)
head(filtconteos)

##################################################################################################################
##################################################################################################################
# EdgeR

## Clase DGEList
d <- DGEList(counts = filtconteos[,1:6], group = targets$tratamiento) ## Normalization
colnames(d) <- targets$label

## Normalization
d <- calcNormFactors(d)
plotMDS(d, main="plotMDS TSvsTC")

## Dispersors stimation
d <- estimateCommonDisp(d,verbose=TRUE)
d <- estimateTagwiseDisp(d)
plotBCV(d, main="plotBCV TSvsTC")

## Test
et <- exactTest(d,pair=c("TS","TC"))
top<- topTags(et, n= Inf)
hist(top$table$FDR, breaks = 100, main = "Hist FDR TSvsTC")
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
write.table(topSig_export, "../../data/DGE/EdgeR_TSvsTC_FDR_0.05.txt", sep="\t", row.names=T)

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
write.table(resSig2_export, "../../data/DGE/DESeq2_TSvsTC_FDR_0.05.txt", sep="\t", row.names=T)

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
write.table(topSig_export, "../../data/DGE/EdgeR_TSvsTC_FDR_5.txt", sep="\t", row.names=T)

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
write.table(resSig2_export, "../../data/DGE/DESeq2_TSvsTC_FDR_5.txt", sep="\t", row.names=T)


###################################################
### How many common DE genes exist edgeR vs DESeq2
###################################################

genesDEcomunes <- intersect(genesDEedgeR,genesDEDESeq2) 
head(genesDEcomunes)
str(genesDEcomunes)


