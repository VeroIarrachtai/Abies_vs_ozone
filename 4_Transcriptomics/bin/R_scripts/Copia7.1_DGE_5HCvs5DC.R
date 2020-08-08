
library(VennDiagram)
library(limma)
library(edgeR)
library(DESeq2)
library(ggbiplot)
library (ggplot2)

# Convert dataframe to data matrix
alldata<-read.delim("../../metadata/all_genes/allreadsgenes.txt")
alldata <- as.data.frame(alldata)
  
x<-alldata
rownames(x)<-alldata[,1]
x<-x[ ,2:ncol(x)]

alldata<-as.matrix(x)
dim (alldata)
########################################### Healthy vs Damaged 170 ppb###########################################
##################################################################################################################
##################################################################################################################
# Select data to compare(descart data)
DCvsHC<- subset(alldata, select = -c(DS_1, DS_2, DS_4,
                                     HS_1, HS_2, HS_5,
                                     HC17, DC47))
alldata<- DCvsHC

############################################
# Add characteristics 
############################################
tratamiento <- c("DC","DC","DC","DC","DC",
                 "HC","HC","HC","HC","HC")
label<- c("DC_1", "DC_2","DC_3","DC_4","DC_5",
          "HC_1","HC_2","HC_3","HC_4","HC_5")
samples<-c("DC1", "DC2","DC3","DC4","DC5",
           "HC1","HC2","HC3","HC4","HC5")
targets<- data.frame(tratamiento,label,samples)
rownames(targets)<- label

targets


# Creamos el objeto con la clase DGEList
y <- DGEList(counts=alldata[,1:10], group=targets$Treatment)

colnames(y) <- targets$Label # Le ponemos nombres a las columnas.
dim(y)
head(y)
### Filtering genes 
keep <- rowSums(cpm(y)>1) >= 5

y <- y[keep,] # Nuevo objeto de estudio, ya filtrado. 
dim(y) # Nueva dimensión

y$samples$lib.size # Tamaño de librería antes del filtrado [1] 978576 1156844 1442169 1485604 1823460 1834335 681743
y$samples$lib.size <- colSums(y$counts)
y$samples$lib.size # Nuevo tamaño de la librería

y <- calcNormFactors(y) 
y$samples

y <- estimateCommonDisp(y,verbose=TRUE) 

y <- estimateTagwiseDisp(y)
plotBCV(y)

et <- exactTest(y) 
et
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
res <- results(dds)

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
topTags(et, n= Inf)

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
genesDEedgeR <- rownames(topSig)
genesDEedgeR
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
write.table(topSig_export, "../../metadata/DGE/EdgeR_HvsD170ppb_FDR_0.05.txt", sep="\t", row.names=T)

topSig <- top[top$table$FDR < 5, ]
dim(topSig)
genesDEedgeR <- rownames(topSig)
genesDEedgeR
topSig_export<-topSig
topSig_export$ID<-genesDEedgeR
write.table(topSig_export, "../../metadata/DGE/EdgeR_HvsD170ppb_FDR_5.txt", sep="\t", row.names=T)


########
# DESeq2
########
# Sort by p-valores
resOrdered <- res[order(res$padj),]
# Only DEG
xx <-res[order(res$padj,na.last=NA),] 
resSig2 <- xx[xx$padj < 0.05, ]
dim(resSig2)
genesDEDESeq2 <- rownames(resSig2)
genesDEDESeq2 
resSig2_export<-resSig2
resSig2_export$ID<-resSig2_export
write.table(resSig2_export, "../../metadata/DGE/DESeq2_HvsD170ppb_FDR_0.05.txt", sep="\t", row.names=T)

xx <-res[order(res$padj,na.last=NA),] 
resSig2 <- xx[xx$padj < 5, ]
dim(resSig2)
genesDEDESeq2 <- rownames(resSig2)
genesDEDESeq2 
resSig2_export<-resSig2
resSig2_export$ID<-resSig2_export
write.table(resSig2_export, "../../metadata/DGE/DESeq2_HvsD170ppb_FDR_5.txt", sep="\t", row.names=T)

###################################################
### How many common DE genes exist edgeR vs DESeq2
###################################################
genesDEcomunes <- intersect(genesDEedgeR,genesDEDESeq2) 
head(genesDEcomunes)
str(genesDEcomunes)

############################
### Plotear diagrama de Venn
############################
grid.newpage() ##To clean the graphics window
dim(resSig2)
dim(topSig)
str(genesDEcomunes)
plot2 <- draw.pairwise.venn(22,20,7,category = c("DESeq2","edgeR"),
                            lty = "blank", 
                            fill = c("cyan3", "hotpink2"))
genesDEcomunes

