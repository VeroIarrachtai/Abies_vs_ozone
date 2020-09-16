# Make Volcano plot
# Verónica Reyes 
# agosto 2020
 
library(limma)
library(edgeR)
library(DESeq2)
library(ggbiplot)
library (ggplot2)


#Load data DESeq2 and EdgeR
results_DESeq2<- read.delim("../../data/DGE/DESeq2_TvsD170ppb_FDR_5.txt")
results_Edge<- read.delim("../../data/DGE/EdgeR_TvsD170ppb_FDR_5.txt")

#Indicate data to plot##
results_DESeq2$sig <- -log10(results_DESeq2$padj) ##Create a column with aditional info of FDR (padj)##
cols_D <- densCols(results_DESeq2$log2FoldChange, results_DESeq2$pvalue)

results_Edge$sig <- -log10(results_Edge$FDR) ##Create a column with aditional info of FDR (padj)##
cols_E <- densCols(results_Edge$logFC, results_Edge$PValue)

##Indicate color code##

cols_D[results_DESeq2$log2FoldChange < -1.5] <- "#0066FF"
cols_D[results_DESeq2$log2FoldChange > 1.5] <- "#0033CC"
cols_D[results_DESeq2$pvalue == 0] <- "#000000"
cols_D[results_DESeq2$pvalue > 0.05] <- "#CCCCCC"
cols_D[results_DESeq2$padj < 0.05] <- "black"
results_DESeq2$cols<- cols_D

cols_E[results_Edge$logFC < -1.5] <- "#0066FF"
cols_E[results_Edge$logFC > 1.5] <- "#0033CC"
cols_E[results_Edge$PValue == 0] <- "#000000"
cols_E[results_Edge$PValue > 0.05] <- "#CCCCCC"
results_Edge$cols<- cols_E

##For the new solution##

VPSol_DESeq2 <- data.frame(results_DESeq2$log2FoldChange, results_DESeq2$sig, row.names = rownames(results_DESeq2)) ##Creates a data frame with coordinates##

VPSol_edge <- data.frame(results_Edge$logFC, results_Edge$sig, row.names = rownames(results_Edge)) ##Creates a data frame with coordinates##


##Rename the columns to simplify the dataframe##

colnames(VPSol_DESeq2) <- c("FoldChange", "pValue")

colnames(VPSol_edge) <- c("FoldChange", "pValue")

##Create a column with colors depending on the value of Fold Change and p-value##  

VPSol_DESeq2$color <- ifelse((VPSol_DESeq2$FoldChange > 1) & (VPSol_DESeq2$pValue < 0.05), "Col_1",
                             ifelse((VPSol_DESeq2$FoldChange < -1) & (VPSol_DESeq2$pValue  < 0.05), "Col_2",
                                    ifelse((VPSol_DESeq2$FoldChange > 1) & (VPSol_DESeq2$pValue  > 0.05), "Col_3",
                                           ifelse((VPSol_DESeq2$FoldChange < -1) & (VPSol_DESeq2$pValue  > 0.05), "Col_4",
                                                  ifelse((VPSol_DESeq2$FoldChange < 1) & (VPSol_DESeq2$pValue  > 0.05), "Col_5", "Col_6")))))

VPSol_DESeq2col_1_sobre <- VPSol_DESeq2[VPSol_DESeq2$color == "Col_1",]
write.table(VPSol_DESeq2col_1_sobre, "../../data/DGE/VPSol_DESeq2col_1_sobre.txt", sep="\t", row.names=T)
VPSol_DESeq2col_2_sub <- VPSol_DESeq2[VPSol_DESeq2$color == "Col_2",]
write.table(VPSol_DESeq2col_2_sub, "../../data/DGE/VPSol_DESeq2col_2_sub.txt", sep="\t", row.names=T)


VPSol_edge$color <- ifelse((VPSol_edge$FoldChange > 1) & (VPSol_edge$pValue  < 0.05), "Col_1",
                           ifelse((VPSol_edge$FoldChange < -1) & (VPSol_edge$pValue  < 0.05), "Col_2",
                                  ifelse((VPSol_edge$FoldChange > 1) & (VPSol_edge$pValue  > 0.05), "Col_3",
                                         ifelse((VPSol_edge$FoldChange < -1) & (VPSol_edge$pValue  > 0.05), "Col_4",
                                                ifelse((VPSol_edge$FoldChange < 1) & (VPSol_edge$pValue  > 0.05), "Col_5", "Col_6")))))

VPSol_edgecol_1_sobre <- VPSol_edge[VPSol_edge$color == "Col_1",]
write.table(VPSol_edgecol_1_sobre, "../../data/DGE/VPSol_edgecol_1_sobre.txt", sep="\t", row.names=T)

VPSol_edgecol_2_sub <- VPSol_edge[VPSol_edge$color == "Col_2",]
write.table(VPSol_edgecol_2_sub, "../../data/DGE/VPSol_edgecol_2_sub.txt", sep="\t", row.names=T)

##Create plot##
ggplot(VPSol_DESeq2, aes(x=FoldChange, y=pValue)) +
  geom_point(aes(colour = color ))
ggsave("../../outputs/8.1_VPSol_DESeq2_FDR0.05.png")

ggplot(VPSol_edge, aes(x=FoldChange, y=pValue)) +
  geom_point(aes(colour = color))
ggsave("../../outputs/8.1_VPSol_edge_FDR0.05.png")





#Load data
results_DESeq2<- read.delim("../../data/DGE/DESeq2_TvsD170ppb_FDR_5.txt")
results_Edge<- read.delim("../../data/DGE/EdgeR_TvsD170ppb_FDR_5.txt")

#Indicate data to plot##
results_DESeq2$sig <- -log10(results_DESeq2$padj) ##Create a column with aditional info of FDR (padj)##
cols <- densCols(results_DESeq2$log2FoldChange, results_DESeq2$pvalue)

results_Edge$sig <- -log10(results_Edge$FDR) ##Create a column with aditional info of FDR (padj)##
cols <- densCols(results_Edge$logFC, results_Edge$PValue)


##Indicate color code##

cols[results_DESeq2$log2FoldChange < -1.5] <- "#0066FF"
cols[results_DESeq2$log2FoldChange > 1.5] <- "#0033CC"
cols[results_DESeq2$pvalue == 0] <- "#000000"
cols[results_DESeq2$pvalue > 0.05] <- "#CCCCCC"


cols[results_Edge$logFC < -1.5] <- "#0066FF"
cols[results_Edge$logFC > 1.5] <- "#0033CC"
cols[results_Edge$PValue == 0] <- "#000000"
cols[results_Edge$PValue > 0.05] <- "#CCCCCC"



##Other graphical parameters##

results_DESeq2$pch <- 19
results_DESeq2$pch[results_DESeq2$pvalue ==0] <- 6
plot(results_DESeq2$log2FoldChange,
     results_DESeq2$sig,
     col=cols, panel.first=grid(),
     main="all RNA",
     xlab="log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=results_DESeq2$pch, cex=0.8)
abline(v=0)
abline(v=c(-1,1), col="brown")


results_Edge$pch <- 19
results_Edge$pch[results_Edge$PValue ==0] <- 6
plot(results_Edge$logFC,
     results_Edge$sig,
     col=cols, panel.first=grid(),
     main="all RNA",
     xlab="log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=results_Edge$pch, cex=0.8)
abline(v=0)
abline(v=c(-1,1), col="brown")


##For the new solution##

VPSol_DESeq2 <- data.frame(results_DESeq2$log2FoldChange, results_DESeq2$sig, row.names = rownames(results_DESeq2)) ##Creates a data frame with coordinates##

VPSol_edge <- data.frame(results_Edge$logFC, results_Edge$sig, row.names = rownames(results_Edge)) ##Creates a data frame with coordinates##


##Rename the columns to simplify the dataframe##

colnames(VPSol_DESeq2) <- c("FoldChange", "pValue")

colnames(VPSol_edge) <- c("FoldChange", "pValue")

##Create a column with colors depending on the value of Fold Change and p-value##  

VPSol_DESeq2$color <- ifelse((VPSol_DESeq2$FoldChange > 1) & (VPSol_DESeq2$pValue < 0.05), "Col_1",
                      ifelse((VPSol_DESeq2$FoldChange < -1) & (VPSol_DESeq2$pValue  < 0.05), "Col_2",
                             ifelse((VPSol_DESeq2$FoldChange > 1) & (VPSol_DESeq2$pValue  > 0.05), "Col_3",
                                    ifelse((VPSol_DESeq2$FoldChange < -1) & (VPSol_DESeq2$pValue  > 0.05), "Col_4",
                                           ifelse((VPSol_DESeq2$FoldChange < 1) & (VPSol_DESeq2$pValue  > 0.05), "Col_5", "Col_6")))))

VPSol_DESeq2col_1_sobre <- VPSol_DESeq2[VPSol_DESeq2$color == "Col_1",]
write.table(VPSol_DESeq2col_1_sobre, "../../metadata/DGE/VPSol_DESeq2col_1_sobre.txt", sep="\t", row.names=T)
VPSol_DESeq2col_2_sub <- VPSol_DESeq2[VPSol_DESeq2$color == "Col_2",]
write.table(VPSol_DESeq2col_2_sub, "../../metadata/DGE/VPSol_DESeq2col_2_sub.txt", sep="\t", row.names=T)


VPSol_edge$color <- ifelse((VPSol_edge$FoldChange > 1) & (VPSol_edge$pValue  < 0.05), "Col_1",
                      ifelse((VPSol_edge$FoldChange < -1) & (VPSol_edge$pValue  < 0.05), "Col_2",
                             ifelse((VPSol_edge$FoldChange > 1) & (VPSol_edge$pValue  > 0.05), "Col_3",
                                    ifelse((VPSol_edge$FoldChange < -1) & (VPSol_edge$pValue  > 0.05), "Col_4",
                                           ifelse((VPSol_edge$FoldChange < 1) & (VPSol_edge$pValue  > 0.05), "Col_5", "Col_6")))))
VPSol_edgecol_1_sobre <- VPSol_edge[VPSol_edge$color == "Col_1",]
write.table(VPSol_edgecol_1_sobre, "../../metadata/DGE/VPSol_edgecol_1_sobre.txt", sep="\t", row.names=T)

VPSol_edgecol_2_sub <- VPSol_edge[VPSol_edge$color == "Col_2",]
write.table(VPSol_edgecol_2_sub, "../../metadata/DGE/VPSol_edgecol_2_sub.txt", sep="\t", row.names=T)

##Create plot##
ggplot(VPSol_DESeq2, aes(x=FoldChange, y=pValue)) +
  geom_point(aes(colour = color ))
ggsave("../../outputs/8.1_VPSol_DESeq2_FDR5.png")

ggplot(VPSol_edge, aes(x=FoldChange, y=pValue)) +
  geom_point(aes(colour = color))
ggsave("../../outputs/8.1_VPSol_edge_FDR5.png")

