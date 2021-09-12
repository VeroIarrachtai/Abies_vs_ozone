# PCA
# Septiembre 2019
# Verónica Reyes


# Load libraries
library (ggplot2)
library("devtools")
library("ggbiplot")

# Load data
#Dar nombre a la base de datos y adjuntarla
metabolites<-read.delim("../../metadata/calculate_relative_abs.txt")

rownames(metabolites)<-c("ST15_1", "ST15_2", "ST15_3", "ST15_4", "ST15_5",
                         "ST16_1", "ST16_2", "ST16_3", "ST16_4", "ST16_5", 
                         "SD15_1", "SD15_2", "SD15_3", "SD15_4", "SD15_5",
                         "SD16_1", "SD16_2", "SD16_3", "SD16_4", "SD16_5",
                         "CS15_1", "CS15_2", "CS15_3", "CS15_4", "CS15_5",
                         "CS16_1", "CS16_2", "CS16_3", "CS16_4", "CS16_5",
                         "CD15_1", "CD15_2", "CD15_3", "CD15_4", "CD15_5",
                         "CD16_1", "CD16_2", "CD16_3", "CD16_4" ,"CD16_5")

beta <- intToUtf8(946)
alfa <-intToUtf8(945)
delta<-intToUtf8(948)
beta.pinene <- paste0(beta,"-Pinene") 
L.alfa.bornyl.acetate <- paste0("L-",alfa,"-Bornyl acetate") 
beta.Caryophyllene.oxide <- paste0(beta,"-Caryophyllene oxide") 
alfa.Caryophyllene <- paste0(alfa,"-Caryophyllene") 
beta.Cubebene <- paste0(beta,"-Cubebene") 
alfa.Cubebene <- paste0(alfa,"-Cubebene") 
delta.Cadinene <- paste0(delta,"-Cadinene") 
alfa.Muurolene <- paste0(alfa,"-Muurolene") 

#moderated period HvsD

metaboliteSTS<-metabolites[1:20,8:14]
colnames(metaboliteSTS)<-c(beta.pinene,L.alfa.bornyl.acetate,beta.Caryophyllene.oxide,
                           alfa.Caryophyllene,beta.Cubebene,alfa.Cubebene,delta.Cadinene)

metabol_ST.pca <- prcomp(metaboliteSTS ,scale.=TRUE)
summary(metabol_ST.pca)
str(metabol_ST.pca)
summary(metabol_ST.pca)
sum<-summary(metabol_ST.pca)
metabolites.PCA<-c(rep("ST15",5), rep("ST16",5), 
                   rep("SD15",5),rep("SD16",5))
metabolites.condition<-c(rep("Tolerant",10), rep("Damaged",10))

ggbiplot(metabol_ST.pca,choices = c(1,2),ellipse=TRUE,obs.scale = 1, var.scale = 1,  groups=metabolites.PCA)+
  scale_color_manual(name="Exposition year",labels = c("2015", "2016", "2015", "2016"), 
                     values=c("#c6003a", "#cdc339", "#7785cc","#59c8a2")) +
  scale_shape_manual(name="Condition", values=c(15,16)) +
  geom_point(aes(colour=metabolites.PCA, shape=metabolites.condition), size = 3)+
  xlab(paste0("Eigenvector 1 explaining ", sum$importance[2,1]*100, "%")) +
  ylab(paste0("Eigenvector 2 explaining ", sum$importance[2,2]*100, "%"))+
  theme_light(base_size = 10)
  
ggsave("../../outputs/5.1_PCA_moderate_TvsD.png")

# "#c6003a", "#cdc339", "#7785cc","#59c8a2"
# "#cd565b", "#e98382", "#00901e","#b1e787"

# #Contingency period HvsD

metabolitesConti<-metabolites[21:40,8:14]
colnames(metabolitesConti)<-c(beta.pinene,L.alfa.bornyl.acetate,beta.Caryophyllene.oxide,
                           alfa.Caryophyllene,beta.Cubebene,alfa.Cubebene,delta.Cadinene)

metabol_Conti.pca <- prcomp(metabolitesConti ,scale.=TRUE)
summary(metabol_Conti.pca)
str(metabol_Conti.pca)
summary(metabol_Conti.pca)
sum<-summary(metabol_Conti.pca)
metabolites.PCA<-c(rep("CS15",5), rep("CS16",5), 
                   rep("CD15",5),rep("CD16",5))
metabolites.condition<-c(rep("Tolerant",10), rep("Damaged",10))

ggbiplot(metabol_Conti.pca,choices = c(1,2),ellipse=TRUE,obs.scale = 1, var.scale = 1,  groups=metabolites.PCA)+
  scale_color_manual(name="Exposition year",labels = c("2015", "2016", "2015", "2016"), 
                     values=c("#c6003a", "#cdc339", "#7785cc","#59c8a2"))+
  scale_shape_manual(name="Condition", values=c(15,16)) +
  geom_point(aes(colour=metabolites.PCA, shape=metabolites.condition), size = 3)+
  xlab(paste0("Eigenvector 1 explaining ", sum$importance[2,1]*100, "%")) +
  ylab(paste0("Eigenvector 2 explaining ", sum$importance[2,2]*100, "%"))+
  theme_light(base_size = 10)
ggsave("../../outputs/5.1_PCA_high_TvsD.png")

# "#c6003a", "#cdc339", "#7785cc","#59c8a2"
# "#cd565b", "#e98382", "#00901e","#b1e787"

####################################################################
# Healthy moderated period vs contingency period

metaboliteshelthy<-metabolites[c(1:10, 21:30),8:14]
colnames(metabolitesConti)<-c(beta.pinene,L.alfa.bornyl.acetate,beta.Caryophyllene.oxide,
                              alfa.Caryophyllene,beta.Cubebene,alfa.Cubebene,delta.Cadinene)

metabolhelthy.pca <- prcomp(metaboliteshelthy ,scale.=TRUE)
summary(metabolhelthy.pca)
str(metabolhelthy.pca)
summary(metabolhelthy.pca)
sum<-summary(metabolhelthy.pca)
metabolites.PCA<-c(rep("ST15",5), rep("ST16",5), 
                   rep("CS15",5),rep("CS16",5))
metabolites.condition<-c(rep("M. concentration",10), rep("H. concentration",10))

ggbiplot(metabolhelthy.pca,choices = c(1,2),ellipse=TRUE,obs.scale = 1, var.scale = 1,  groups=metabolites.PCA)+
  scale_color_manual(name="Exposition year",labels = c("2015", "2016", "2015", "2016"), 
                     values=c("#b35bb7", "#a3bb4b", "#7785cc","#59c8a2"))+
  scale_shape_manual(name="Period", values=c(15,16)) +
  geom_point(aes(colour=metabolites.PCA, shape=metabolites.condition), size = 3)+
  xlab(paste0("Eigenvector 1 explaining ", sum$importance[2,1]*100, "%")) +
  ylab(paste0("Eigenvector 2 explaining ", sum$importance[2,2]*100, "%"))+
  theme_light(base_size = 10)
ggsave("../../outputs/5.1_PCA_tolerant_modevshigh.png")

# "#b35bb7", "#a3bb4b", "#7785cc","#59c8a2"
# "#6cb643","#b4a945","#53b483","#617835"

#Damaged moderated period vs contingency period

metabolitesdamaged<-metabolites[c(11:20, 31:40),8:14]
colnames(metabolitesdamaged)<-c(beta.pinene,L.alfa.bornyl.acetate,beta.Caryophyllene.oxide,
                              alfa.Caryophyllene,beta.Cubebene,alfa.Cubebene,delta.Cadinene)
metaboldamaged.pca <- prcomp(metabolitesdamaged ,scale.=TRUE)
summary(metaboldamaged.pca)
str(metaboldamaged.pca)
summary(metaboldamaged.pca)
sum<-summary(metaboldamaged.pca)
metabolites.PCA<-c(rep("SD15",5), rep("SD16",5), 
                   rep("CD15",5),rep("CD16",5))
metabolites.condition<-c(rep("M. concentration",10), rep("H. concentration",10))

ggbiplot(metaboldamaged.pca,choices = c(1,2),ellipse=TRUE,obs.scale = 1, var.scale = 1,  groups=metabolites.PCA)+
  scale_color_manual(name="Exposition year",labels = c("2015", "2016", "2015", "2016"), 
                     values=c("#c6003a", "#cdc339", "#8cd743", "#d59832")) +
  scale_shape_manual(name="Period", values=c(15,16)) +
  geom_point(aes(colour=metabolites.PCA, shape=metabolites.condition), size = 3)+
  xlab(paste0("Eigenvector 1 explaining ", sum$importance[2,1]*100, "%")) +
  ylab(paste0("Eigenvector 2 explaining ", sum$importance[2,2]*100, "%"))+
  theme_light(base_size = 10)

ggsave("../../outputs/5.1_PCA_damaged_modevshigh.png")

# "#c6003a", "#cdc339", "#df355f", "#8cd743"
# "#dd5035","#d7ac43","#c2455e","#a95f2e"
