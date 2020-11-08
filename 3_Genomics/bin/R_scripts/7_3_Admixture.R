# Reyes Galindo Verónica
# 12 Febrero 2019

### This script 
# Script proporcionado por Lewis Spruning 

## Load packages
library(magrittr)
library(reshape)
library(ggplot2)

## Load data cross validation
log <- read.table("../../metadata/admixture_PQ_files/logall_snp_withoutDupLoci_88s_maxmiss0.9_maf0.05.txt")
ggplot(data=log, aes(x= V1, y= V2))+
  geom_point(colour = "#ce46ac") +
  geom_line(colour = "#ce46ac") +  
  ylim(.40, 1.2) +
  labs(x="K Value" , y="Cross-validation error")+
  theme(axis.title.x = element_text(size = rel(1.3)))+
  theme(axis.title.y = element_text(size = rel(1.3)))+
  theme(axis.text.x = element_text(hjust = .5, size=11, color="black"))+
  theme(axis.text.y = element_text(hjust = .5, size=11, color="black"))+
  theme_light()

ggsave("../../outputs/7.3_Cross_val.png")

structureplot <- function(str_out,pops,k,xlab = F)
{

#Sort columns by prevalence
str_out <- str_out[,order(apply(str_out,2,sum),decreasing = T)]
  
  x <- melt(as.matrix(str_out)) %>%
    transform(X1 = factor(X1, rownames(str_out)),
              X2 = factor(X2, colnames(str_out)),
              value = as.numeric(paste(value)))
  
  x$pop <- rep(pops[,1],k)
  str_out$pop <- pops[,1]
  x$X3 <- factor(paste(letters[as.numeric(x$pop)],x$X1,sep = "_"))
  x <- x[order(x$X3),]
  
  pop_pos <- cumsum(tapply(1:nrow(str_out),str_out$pop,length))
  labpos <- pop_pos
  
  for(i in 1:length(pop_pos))
  {
    if(i == 1)
    {
      labpos[i] <- pop_pos[i]/2
    } else
    {
      labpos[i] <- pop_pos[i-1] + (pop_pos[i]-pop_pos[i-1])/2
    }
  }
  labpos <- pop_pos
  
  for(i in 1:length(pop_pos))
  {
    if(i == 1)
    {
      labpos[i] <- pop_pos[i]/2
    } else
    {
      labpos[i] <- pop_pos[i-1] + (pop_pos[i]-pop_pos[i-1])/2
    }
  }
  
  labpos <- round(labpos,0)
  
  #Specify some nice colours
  mycols <- c("#c167ba",
              "#93c251",
              "#716cd0",
              "#cc9933",
              "#c75b81",
              "#5dc397"
  )[1:k]
  
  if(xlab ==T)
  {
    output <- ggplot(x,
                     aes(x = as.numeric(X3),y = value, fill = X2))+
      geom_bar(stat = "identity",width = 2)+
      theme_bw()+
      scale_x_continuous(breaks = labpos,labels = names(pop_pos),expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      geom_vline(xintercept = pop_pos)+
      xlab("")+
      ylab("")+
      theme(axis.text.x = element_text(vjust = 0.1,hjust = 1,angle = 90),
            axis.ticks = element_blank(),
            legend.position = "none")+
      scale_fill_manual(values = mycols)
  } else
  {
    output <- ggplot(x,
                     aes(x = as.numeric(X3),y = value, fill = X2))+
      geom_bar(stat = "identity",width = 2)+
      theme_bw()+
      scale_x_continuous(breaks = labpos,labels = rep("",length(pop_pos)),expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      geom_vline(xintercept = pop_pos)+
      xlab("")+
      ylab("")+
      theme(axis.text.x = element_text(vjust = 0.1,hjust = 1,angle = 90),
            axis.ticks = element_blank(),
            legend.position = "none")+
      scale_fill_manual(values = mycols)
  }
  
  
  return(output)
}


structureplot_popord <- function(str_out,pops,k,xlab = T, target_order)
{
  #Sort columns by prevalence
  str_out <- str_out[,order(apply(str_out,2,sum),decreasing = T)]
  
  x <- melt(as.matrix(str_out)) %>%
    transform(X1 = factor(X1, rownames(str_out)),
              X2 = factor(X2, colnames(str_out)),
              value = as.numeric(paste(value)))
  
  x$pop <- rep(pops[,1],k)
  str_out$pop <- pops[,1]
  
  #reorder x by pop in target order
  
  x<-x[order(unlist(sapply(x$pop, function(y) which(target_order == y)))),]
  x$pop<-factor(x$pop, levels=target_order)
  str_out$pop <-factor(str_out$pop, levels=target_order)
  
  # add pop codes and reorder
  x$X3 <- factor(paste(letters[as.numeric(x$pop)],x$X1,sep = "_"))
  x <- x[order(x$X3),]
  
  pop_pos <- cumsum(tapply(1:nrow(str_out),str_out$pop,length))
  labpos <- pop_pos
  
  for(i in 1:length(pop_pos))
  {
    if(i == 1)
    {
      labpos[i] <- pop_pos[i]/2
    } else
    {
      labpos[i] <- pop_pos[i-1] + (pop_pos[i]-pop_pos[i-1])/2
    }
  }
  
  labpos <- round(labpos,0)
  
  
  #Specify some nice colours
  mycols <- c("#c167ba",
              "#93c251",
              "#716cd0",
              "#cc9933",
              "#c75b81",
              "#5dc397"
  )[1:k]
  
  if(xlab ==T)
  {
    output <- ggplot(x,
                     aes(x = as.numeric(X3),y = value, fill = X2))+
      geom_bar(stat = "identity",width = 2)+
      theme_bw()+
      scale_x_continuous(breaks = labpos,labels = names(pop_pos),expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      geom_vline(xintercept = pop_pos)+
      xlab("")+
      ylab("")+
      theme(axis.text.x = element_text(vjust = 0.1,hjust = 1,angle = 90),
            axis.ticks = element_blank(),
            legend.position = "none")+
      scale_fill_manual(values = mycols)
  } else
  {
    output <- ggplot(x,
                     aes(x = as.numeric(X3),y = value, fill = X2))+
      geom_bar(stat = "identity",width = 2)+
      theme_bw()+
      scale_x_continuous(breaks = labpos,labels = rep("",length(pop_pos)),expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      geom_vline(xintercept = pop_pos)+
      xlab("")+
      ylab("")+
      theme(axis.text.x = element_text(vjust = 0.1,hjust = 1,angle = 90),
            axis.ticks = element_blank(),
            legend.position = "none")+
      scale_fill_manual(values = mycols)
  }
  
  
  return(output)
}



## Load data about family
famfile <- read.table("../../data/without_Dup_loci/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05.fam")
par(mfrow=c(4,1))


## Draw differents K´s MODIFIED FUN ALICIA

z<-c("SierraManantlan","NevadoColima","VolcanTancitaro","PuertaGarnica","MichoacanAlt",
     "SanAndres","CerroZamorano","CerroBlanco","NevadoTolucaRG","NevadoTolucaNT",
     "NevadoTolucaSB","SantaRosaXochiacD", "SantaRosaXochiacS", "Ajusco","Ixtapalucan",
     "ElChico","Tlaxco","VolcanAtlitzin","Malinche","CofrePerote")

# K2
qfile <- read.table("../../metadata/admixture_PQ_files/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05.2.Q")
structureplot_popord(str_out = qfile, pops = famfile, k = 2, xlab = T, target_order=z)
ggsave("../../outputs/7.3_Admixture_K2.png")

# K3
qfile <- read.table("../../metadata/admixture_PQ_files/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05.3.Q")
structureplot_popord(str_out = qfile,pops = famfile,k = 3, xlab = T, target_order =z)
ggsave("../../outputs/7.3_Admixture_K3.png")

# K4
qfile <- read.table("../../metadata/admixture_PQ_files/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05.4.Q")
structureplot_popord(str_out = qfile,pops = famfile,k = 4, xlab = T, target_order =z)
ggsave("../../outputs/7.3_Admixture_K4.png")

# K5
qfile <- read.table("../../metadata/admixture_PQ_files/snp_withoutDupLoci_88s_maxmiss0.9_maf0.05.5.Q")
structureplot_popord(str_out = qfile,pops = famfile,k = 5, xlab = T, target_order =z)
ggsave("../../outputs/7.3_Admixture_K5.png")
