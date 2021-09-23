# Verónica Reyes Galindo´
# febrero 2020
# merge figures

# Load data

C_TvsD<-read.delim("../../data/Over_Down/GENERAL_C_TvsD.txt", sep="\t")
D_170Cvs87SS<-read.delim("../../data/Over_Down/GENERAL_D_170Cvs87SS.txt", sep="\t")
T_170Cvs87SS<-read.delim("../../data/Over_Down/GENERAL_T_170Cvs87SS.txt", sep="\t")

# Plot individual

ggplot(C_TvsD, aes(x=log2FoldChange_D2, y=sig_D2)) +
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

ggplot(D_170Cvs87SS, aes(x=log2FoldChange_D2, y=sig_D2)) +
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


ggplot(T_170Cvs87SS, aes(x=log2FoldChange_D2, y=sig_D2)) +
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

# Merge plots

library(ggplot2)
library(ggpubr)


