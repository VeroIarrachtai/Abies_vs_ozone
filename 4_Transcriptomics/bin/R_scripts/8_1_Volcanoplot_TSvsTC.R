# Verónica Reyes Galindo

#Load libraries
library(ggplot2)
library(ggrepel)

#Load data DESeq2 and EdgeR
results_DESeq2<- read.delim("../../data/DGE/DESeq2_TSvsTC_FDR_5.txt") ## TENGO ESTO 2 VECES EN LA MISMA TABLA
results_DESeq2<- results_DESeq2[,1:6]
results_Edge<- read.delim("../../data/DGE/EdgeR_TSvsTC_FDR_5.txt")

# Change value EdgeR LoadChange

results_Edge$logFC <- -(results_Edge$logFC)


# Remake tables with info that I want
results_DESeq2_rt <- data.frame(results_DESeq2[,c(2,5,6)])
results_Edge_rt <- data.frame(results_Edge[,c(1,3,4)])

results_DESeq2_rt$sig <- -log10(results_DESeq2$pvalue) ##Create a column with aditional info of FDR (padj)##
results_DESeq2_rt$sigadj <- -log10(results_DESeq2$padj) ##Create a column with aditional info of FDR (padj)##
results_Edge_rt$sig <- -log10(results_Edge$PValue) ##Create a column with aditional info of FDR (padj)##
results_Edge_rt$sigadj <- -log10(results_Edge$FDR) ##Create a column with aditional info of FDR (padj)##

##Create a column with colors depending p-adj##  
results_DESeq2_rt$TDE <- (results_DESeq2_rt$padj< 0.05) & (abs(results_DESeq2_rt$log2FoldChange) > 1)

results_Edge_rt$TDE <- (results_Edge_rt$FDR< 0.05) & (abs(results_Edge_rt$logFC) > 1)

# Plot Volcano plot
ggplot(results_DESeq2_rt, aes(x=log2FoldChange, y=sig)) +
  geom_point(aes(colour = TDE ),size =3.5)+
  geom_text(aes(label=ifelse((padj< 0.05) & (abs(log2FoldChange) > 1),
                             as.character(row.names(results_DESeq2_rt)),'')),hjust=0.5,vjust=0.5, size= 2, angle=35)+
  scale_color_manual(values=c("grey", "#d44792"))+
  xlab("log2 fold change")+
  ylab("-log10 (P value)")+
  theme_light(base_size = 10)+
  theme(legend.position = "none")+
  geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed", size = 0.25) + # Horizontal significance cut-off line.
  geom_vline(xintercept = 1, colour = "black", linetype = "dashed", size = 0.25)+  # Vertical significance cut-off line (+).
  geom_vline(xintercept = -1, colour = "black", linetype = "dashed", size = 0.25)  # Vertical significance cut-off line (+).
ggsave("../../outputs/8_1_VP_DE2_TSvsTC.png")

ggplot(results_Edge_rt, aes(x=logFC, y=sig)) +
  geom_point(aes(colour = TDE ),size =3.5)+
  geom_text(aes(label=ifelse((FDR< 0.05) & (abs(logFC) > 1),
                             as.character(row.names(results_Edge_rt)),'')),hjust=0.5,vjust=0.5, size= 2, angle=35)+
  scale_color_manual(values=c("grey", "#ddaee8"))+
  xlab("log2 fold change")+
  ylab("-log10 (P value)")+
  theme_light(base_size = 10)+
  theme(legend.position = "none")+
  geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed", size = 0.25) + # Horizontal significance cut-off line.
  geom_vline(xintercept = 1, colour = "black", linetype = "dashed", size = 0.25)+  # Vertical significance cut-off line (+).
  geom_vline(xintercept = -1, colour = "black", linetype = "dashed", size = 0.25)  # Vertical significance cut-off line (+).
ggsave("../../outputs/8_1_VP_ER_TSvsTC.png")

# Make table with over and down expressed
overxpress_Ds2 <- results_DESeq2_rt[results_DESeq2_rt[,"padj"]< 0.05 & results_DESeq2_rt[,"log2FoldChange"]> 1, ]
downxpress_Ds2 <- results_DESeq2_rt[results_DESeq2_rt[,"padj"]< 0.05 & results_DESeq2_rt[,"log2FoldChange"]< -1, ]

overxpress_Ed <- results_Edge_rt[results_Edge_rt[,"FDR"]< 0.05 & results_Edge_rt[,"logFC"]> 1, ]
downxpress_Ed <- results_Edge_rt[results_Edge_rt[,"FDR"]< 0.05 & results_Edge_rt[,"logFC"]< -1, ]

# Obtain comun genes

genesDEcomun_over <- intersect(rownames(overxpress_Ds2),rownames(overxpress_Ed))
genesDEcomun_down <- intersect(rownames(downxpress_Ds2),rownames(downxpress_Ed))

comun_genes_od <- c(genesDEcomun_over,genesDEcomun_down)

# Export data 
write.table(overxpress_Ds2, "../../data/Over_Down/over_DE2_TSvsTC.txt", sep="\t", row.names=T)
write.table(downxpress_Ds2, "../../data/Over_Down/down_DE2_TSvsTC.txt", sep="\t", row.names=T)

write.table(overxpress_Ed , "../../data/Over_Down/over_ER_TSvsTC.txt", sep="\t", row.names=T)
write.table(downxpress_Ed , "../../data/Over_Down/down_ER_TSvsTC.txt", sep="\t", row.names=T)

