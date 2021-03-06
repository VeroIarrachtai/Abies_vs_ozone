# Verónica Reyes Galindo
# febrero 2020
# Damaged vs Tolerant 170 ppb

#Load libraries
library(ggplot2)
library(ggrepel)

#Load data DESeq2 and EdgeR
results_DESeq2<- read.delim("../../data/DGE/DESeq2_TvsD170ppb_FDR_5.txt") ## TENGO ESTO 2 VECES EN LA MISMA TABLA
results_DESeq2<- results_DESeq2[,1:6]
results_Edge<- read.delim("../../data/DGE/EdgeR_TvsD170ppb_FDR_5.txt")

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


# Make table with over and down expressed in TOLERANTS
overxpress_Ds2 <- results_DESeq2_rt[results_DESeq2_rt[,"padj"]< 0.05 & results_DESeq2_rt[,"log2FoldChange"]> 1, ]
downxpress_Ds2 <- results_DESeq2_rt[results_DESeq2_rt[,"padj"]< 0.05 & results_DESeq2_rt[,"log2FoldChange"]< -1, ]

overxpress_Ed <- results_Edge_rt[results_Edge_rt[,"FDR"]< 0.05 & results_Edge_rt[,"logFC"]> 1, ]
downxpress_Ed <- results_Edge_rt[results_Edge_rt[,"FDR"]< 0.05 & results_Edge_rt[,"logFC"]< -1, ]

# Obtain comun genes
genesDEcomun_over <- intersect(rownames(overxpress_Ds2),rownames(overxpress_Ed))
genesDEcomun_down <- intersect(rownames(downxpress_Ds2),rownames(downxpress_Ed))

comun_genes_od <- c(genesDEcomun_over,genesDEcomun_down)
comun_genes_od

# PLOT GENERAL
# To create column with row name
results_DESeq2_rt$rownames <- rownames(results_DESeq2_rt) 
results_Edge_rt$rownames <- rownames(results_Edge_rt)

# To create new data frame with differential expression data (DESeq2 and edgeR) 
df_general<- merge(results_DESeq2_rt, results_Edge_rt, by = "rownames",  all.x=TRUE)
colnames(df_general)<- c( "rownames","log2FoldChange_D2","pvalue_D2","padj_D2" ,
                          "sig_D2", "sigadj_D2","TDE_D2","logFC_ER","PValue_ER",
                          "FDR_ER", "sig_ER", "sigadj_ER","TDE_ER")

###Colors
#Add column with color cathegory accord to diferential expression with DESeq2 and edgeR

df_general$color <- ifelse((df_general$TDE_D2 == TRUE) & (df_general$TDE_ER == TRUE), "Col_1", #DESeq2 and EdgeR
                   ifelse((df_general$TDE_D2 == TRUE) & (df_general$TDE_ER == FALSE), "Col_2", # Only DESeq2 
                          ifelse((df_general$TDE_D2 == FALSE) & (df_general$TDE_ER == TRUE), "Col_3", # Only EdgeR
                                 ifelse((df_general$TDE_D2 == FALSE) & (df_general$TDE_ER == FALSE), "Col_4", "Col_5"))))# No Differential Expression
rownames(df_general)<- df_general$rownames

# Export data
write.table(df_general, "../../data/Over_Down/GENERAL_C_TvsD.txt", sep="\t", row.names=T)


# Plot Volcano plot
ggplot(df_general, aes(x=log2FoldChange_D2, y=sig_D2)) +
  geom_point(aes(colour =  color ),size =3.5)+
  scale_color_manual(values=c("#cd77ca", # purple D2 and ER
                              "#49cfd8", # blue Only D2
                              "#cdbd3a", # yellow Only ER
                              "grey",    # grey
                              "black"))+ # black
  xlab("log2 fold change")+
  ylab("-log10 (P value)")+
  theme_light(base_size = 10)+
  theme(legend.position = "none")+
  geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed", size = 0.25) + # Horizontal significance cut-off line.
  geom_vline(xintercept = 1, colour = "black", linetype = "dashed", size = 0.25)+  # Vertical significance cut-off line (+).
  geom_vline(xintercept = -1, colour = "black", linetype = "dashed", size = 0.25)  # Vertical significance cut-off line (+).
  ggsave("../../outputs/8_1_VP_General_sinN_C_TvsD.png")

# Export data 

# color 1= pink D2 and ER
# color 2= blue Only D2
# color 3= purple Only ER

D2_ER_genes <-subset(df_general, subset= color == "Col_1")
D2_ER_genes$rownames
D2_genes <-subset(df_general, subset= color == "Col_2")
D2_genes$rownames
ER_genes <-subset(df_general, subset= color == "Col_3")
ER_genes$rownames

write.table(D2_ER_genes, "../../data/Over_Down/SPECIFIC/D2_ER_C_TvsD.txt", sep="\t", row.names=T)
write.table(D2_genes, "../../data/Over_Down/SPECIFIC/D2_C_TvsD.txt", sep="\t", row.names=T)
write.table(ER_genes, "../../data/Over_Down/SPECIFIC/ER_C_TvsD.txt", sep="\t", row.names=T)


write.table(genesDEcomun_over ,"../../data/Over_Down/SPECIFIC/IDs_D2_ER_C_TvsD_over.txt",sep = "\t", row.names = F, col.names = F)
write.table(genesDEcomun_down ,"../../data/Over_Down/SPECIFIC/IDs_D2_ER_C_TvsD_down.txt",sep = "\t", row.names = F, col.names = F)

