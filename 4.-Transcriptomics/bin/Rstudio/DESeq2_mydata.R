# Load DESeq2 library
library("DESeq2")


# Load data
alldata<-read.delim("../../allreadsgenes.txt")
alldata <- as.data.frame(alldata)
head (alldata)

# Convert dataframe to data matrix
x<-alldata
rownames(x)<-alldata[,1]
x<-x[ ,2:ncol(x)]
alldata<-as.matrix(x)

########################################### Healthy vs Damaged 170 ppb###########################################
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
label<- c("DC_1", "DC_2","DC_3","DC_4","DC_5",
          "HC_1","HC_2","HC_3","HC_4","HC_5")
samples<-c("DC1", "DC2","DC3","DC4","DC5",
           "HC1","HC2","HC3","HC4","HC5")
targets<- data.frame(tratamiento,label,samples)
rownames(targets)<- label

targets

### Filtering genes 
table(rowSums(DCvsHC)==0)
suma <- rowSums(DCvsHC)
filtconteos <- DCvsHC[suma != 0,] 
dim(filtconteos)
head(filtconteos)
# DESeq2

### Class DESeqDataSet
ddsHTSeq <- DESeqDataSetFromMatrix(countData=filtconteos, colData= targets, design=~tratamiento)

# Differential expression analysis
#differential expression analysis steps are wrapped into a single function, DESeq()
dds <- DESeq(ddsHTSeq)
# restuls talbe will be generated using results() which will include:
#  log2 fold changes, p values and adjusted p values
res <- results(dds)
res
summary(res)

resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds,normalized =TRUE)), 
                 by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'

outputPrefix <- "Abies_DESeq2"

write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized_all.csv"))

# filter results by p value
res= subset(res, padj<0.1)

# order results by padj value (most significant to least)
res <- res[order(res$padj),]
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds,normalized =TRUE)), 
                 by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'

outputPrefix <- "Abies_DESeq2"

write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), 
            file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),
          file = paste0(outputPrefix, "-test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
temp_ddsClean <- ddsClean
tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)
resClean = subset(res, padj<0.1)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) PCA plot
# 4) heatmap of clustering analysis
#
#
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off() 

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld),file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd),file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t'))


# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
sampleDists <- dist(t(assay(rld)))
suppressMessages(library("RColorBrewer"))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(rld), rld$type, sep="")
colnames(sampleDistMatrix) <- paste(colnames(rld), rld$type, sep="")
colors <- colorRampPalette( rev(brewer.pal(8, "Blues")) )(255)
heatmap(sampleDistMatrix,col=colors,margin = c(8,8))
dev.copy(png,paste0(outputPrefix, "-clustering.png"))
dev.off()

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
condition <- tratamiento
scores <- data.frame(pc$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))
#dev.copy(png,paste0(outputPrefix, "-PCA.png"))
ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.png"))


# heatmap of data
library("RColorBrewer")
library("gplots")
# 35 top expressed genes with heatmap.2
select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:35]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=35)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Heatmap of 35 DE Genes in ozone Comparison")
dev.copy(png, paste0(outputPrefix, "-HEATMAP.png"))
dev.off()
DEgenes_Edge$rowname
