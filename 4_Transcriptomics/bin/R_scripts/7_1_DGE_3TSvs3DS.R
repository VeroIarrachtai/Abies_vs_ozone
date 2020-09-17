# Verónica Reyes,
# febrero 2020
# Damaged vs Tolerant 87 ppb

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
DSvsTS<- subset(alldata, select = c(DS_1, DS_2, DS_4,
                                     TS_1, TS_2, TS_5))


############################################
# Add characteristics 
############################################
tratamiento <- c("DS","DS","DS",
                 "TS","TS","TS")
label <- c("DS_1", "DS_2","DS_4",
          "TS_1","TS_2","TS_5")
samples <-c("DS1", "DS2","DS4",
            "TS1","TS2","TS5")
targets <- data.frame(tratamiento,label,samples)
rownames(targets) <- label

targets

### Filtering genes 
table(rowSums(DSvsTS)==0)
suma <- rowSums(DSvsTS)
filtconteos <- DSvsTS[suma != 0,] 
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
plotMDS(d, main="plotMDS DSvsTS")

## Dispersors stimation
d <- estimateCommonDisp(d,verbose=TRUE)
d <- estimateTagwiseDisp(d)
plotBCV(d, main="plotBCV DSvsTS")

## Test
et <- exactTest(d,pair=c("DS","TS"))
top<- topTags(et, n= Inf)
hist(top$table$FDR, breaks = 100, main = "Hist FDR DSvsTS")
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
# EdgeR
########

topSig <- top[top$table$FDR < 0.05, ]
dim(topSig)
genesDEedgeR <- rownames(topSig)
head(genesDEedgeR)
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
head(topSig_export)
write.table(topSig_export, "../../data/DGE/EdgeR_DSvsTS_FDR_0.05.txt", sep="\t", row.names=T)

########
# DESeq2
########
xx <-res[order(res$padj,na.last=NA),] 
resSig2 <- xx[xx$padj < 0.05, ]
dim(resSig2)
head(resSig2)
genesDEDESeq2 <- rownames(resSig2)
head(resSig2)
resSig2_export<-resSig2
resSig2_export$ID<-resSig2_export
head(resSig2_export)
write.table(resSig2_export, "../../data/DGE/DESeq2_DSvsTS_FDR_0.05.txt", sep="\t", row.names=T)


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
# EdgeR
########

topSig <- top[top$table$FDR < 5, ]
dim(topSig)
genesDEedgeR <- rownames(topSig)
head(genesDEedgeR)
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
head(topSig_export)
write.table(topSig_export, "../../data/DGE/EdgeR_DSvsTS_FDR_5.txt", sep="\t", row.names=T)

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
write.table(resSig2_export, "../../data/DGE/DESeq2_DSvsTS_FDR_5.txt", sep="\t", row.names=T)


###################################################
### How many common DE genes exist edgeR vs DESeq2
###################################################

genesDEcomunes <- intersect(genesDEedgeR,genesDEDESeq2) 
head(genesDEcomunes)
str(genesDEcomunes)


