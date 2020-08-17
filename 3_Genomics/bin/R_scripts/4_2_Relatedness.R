# Reyes Galindo Verónica
# 12 Febrero 2019

library(ggplot2)

# Cargar matriz de datos del calculo de relatedness
pair_rel <- read.table("../../data/relatedness/relsnp_snp_withoutDupLoci_89ind_maxmiss0.9_maf0.05.rel", as.is = T)

#Crear una tabla con el nombre de las muestras pareadas y cargar ID de las muestras
ids <- read.table("../../data/relatedness/relsnp_snp_withoutDupLoci_89ind_maxmiss0.9_maf0.05.rel.id", as.is = T)
colnames(pair_rel) <- ids$V2
rownames(pair_rel) <- ids$V2

x <- data.frame(t(combn(ids$V2,2)))
colnames(x) <- c("IND1", "IND2")
x$rel <- NA
x$POP1 <- NA
x$POP2 <- NA

for(i in 1:nrow(x))
  {
  x$rel[i] <- pair_rel[paste(x$IND1[i]), paste(x$IND2[i])]
  x$POP1[i] <- subset(ids,V2 == x$IND1[i])$V1
  x$POP2[i] <- subset(ids,V2 == x$IND2[i])$V1
}

x2 <- subset(x, POP1 == POP2)


x2$facet = factor(x2$POP1, levels = c("SierraManantlan","NevadoColima", "VolcanTancitaro","PuertaGarnica","MichoacanAlt","SanAndres" ,
                                      "CerroZamorano","CerroBlanco", "NevadoTolucaRG","NevadoTolucaNT","NevadoTolucaSB",
                                      "SantaRosaXochiacD","SantaRosaXochiacS","Ajusco","Ixtapalucan","ElChico",
                                      "Tlaxco","VolcanAtlitzin","Malinche","CofrePerote"))

levels(x2$facet) <- c("Sierra Manantlán","Nevado de Colima", "Volcán Tancítaro","Puerta Garnica","Michoacan Alt","San Andrés" ,
                      "Cerro Zamorano","Cerro Blanco", "Nevado Toluca RG","Nevado Toluca NT","Nevado Toluca SB",
                      "Santa Rosa Xochiac Damaged","Santa Rosa Xochiac Tolerant","Ajusco","Ixtapalucan","El Chico",
                      "Tlaxco","Volcán Atlitzin","Malinche","Cofre de Perote")

x2$colplot <- c(rep("Other",9),rep("Damaged",10), rep("Tolerant",10),rep("Other",130))

ggplot(x2,aes(x=rel, color=colplot, fill=colplot))+
  geom_histogram(bins=25 )+
  facet_wrap(~ facet)+
  scale_color_manual(values=c("#c6003a", "#999999","#00901e"))+
  scale_fill_manual(values=c("#c6003a", "#999999","#00901e"))+
  xlab("Relatedness")+
  ylab("Pairwaise count")+
  theme(axis.title.y = element_text(size = rel(2.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(2.5), angle = 360))+
  theme(axis.text.x = element_text(hjust = 0.5, size=8, color="black"))+
  theme(axis.text.y = element_text(hjust = 0.5, size=8, color="black"))+
  theme_light()+
  theme(legend.position = "none")

ggsave("../../outputs/4.2_Relatedness.png")
