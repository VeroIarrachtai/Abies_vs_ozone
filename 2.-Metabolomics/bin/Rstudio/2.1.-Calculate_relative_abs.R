### This script make a dataframe with all samples of GC-SM analisys###
## February 2020
## Verónica Reyes Galindo


# Load data
metabolites<- read.csv("../../metadata/metabolitos_Tesis_Vero.csv")
metabolites2<- data.frame(metabolites$C.A.T,metabolites$C.A, metabolites$Season,metabolites$Condition, metabolites$A.expo, metabolites$Sample, metabolites$mmet)

colnames(metabolites2)<- c("C.A.T", "C.A", "Season",
                           "Condition", "A.expo", "Sample", "mmet") 

# Calcular X1
metabolites2$x1<-(metabolites$microlitres.standar*metabolites$miligrams/metabolites$microlitres)*2/1000


#Calcular X2,X3,X4,X5 para cada metabolito

# beta.pinene
metabolites2$X2_beta.pinene <- metabolites2$x1*metabolites$beta.pinene/metabolites$Area.standar.porc 
metabolites2$X3_beta.pinene <- metabolites2$X2_beta.pinene*metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_beta.pinene <- metabolites2$X3_beta.pinene*100/metabolites$We.needle.grams
metabolites2$X5_beta.pinene <- metabolites2$X4_beta.pinene*1000
# L.alfa.bornyl.acetate
metabolites2$X2_L.alfa.bornyl.acetate <- metabolites2$x1*metabolites$L.alfa.bornyl.acetate/metabolites$Area.standar.porc
metabolites2$X3_L.alfa.bornyl.acetate <- metabolites2$X2_L.alfa.bornyl.acetate*metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_L.alfa.bornyl.acetate <- metabolites2$X3_L.alfa.bornyl.acetate*100/metabolites$We.needle.grams
metabolites2$X5_L.alfa.bornyl.acetate <- metabolites2$X4_L.alfa.bornyl.acetate*1000
# beta.Caryophyllene.oxide
metabolites2$X2_beta.Caryophyllene.oxide <- metabolites2$x1*metabolites$beta.Caryophyllene.oxide/metabolites$Area.standar.porc
metabolites2$X3_beta.Caryophyllene.oxide <- metabolites2$X2_beta.Caryophyllene.oxide*metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_beta.Caryophyllene.oxide <- metabolites2$X3_beta.Caryophyllene.oxide*100/metabolites$We.needle.grams
metabolites2$X5_beta.Caryophyllene.oxide <- metabolites2$X4_beta.Caryophyllene.oxide*1000
#  alfa.Caryophyllene
metabolites2$X2_alfa.Caryophyllene <- metabolites2$x1*metabolites$alfa.Caryophyllene/metabolites$Area.standar.porc
metabolites2$X3_alfa.Caryophyllene <- metabolites2$X2_alfa.Caryophyllene*metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_alfa.Caryophyllene <- metabolites2$X3_alfa.Caryophyllene*100/metabolites$We.needle.grams
metabolites2$X5_alfa.Caryophyllene <- metabolites2$X4_alfa.Caryophyllene*1000
# beta.Cubebene
metabolites2$X2_beta.Cubebene <- metabolites2$x1*metabolites$beta.Cubebene/metabolites$Area.standar.porc
metabolites2$X3_beta.Cubebene <- metabolites2$X2_beta.Cubebene*metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_beta.Cubebene <- metabolites2$X3_beta.Cubebene*100/metabolites$We.needle.grams
metabolites2$X5_beta.Cubebene <- metabolites2$X4_beta.Cubebene*1000
# alfa.Cubenene
metabolites2$X2_alfa.Cubenene  <- metabolites2$x1*metabolites$alfa.Cubenene /metabolites$Area.standar.porc
metabolites2$X3_alfa.Cubenene  <- metabolites2$X2_alfa.Cubenene *metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_alfa.Cubenene  <- metabolites2$X3_alfa.Cubenene *100/metabolites$We.needle.grams
metabolites2$X5_alfa.Cubenene  <- metabolites2$X4_alfa.Cubenene *1000
# delta.Cadinene
metabolites2$X2_delta.Cadinene <- metabolites2$x1*metabolites$delta.Cadinene/metabolites$Area.standar.porc
metabolites2$X3_delta.Cadinene <- metabolites2$X2_delta.Cadinene*metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_delta.Cadinene <- metabolites2$X3_delta.Cadinene*100/metabolites$We.needle.grams
metabolites2$X5_delta.Cadinene <- metabolites2$X4_delta.Cadinene*1000
# alfa.Muurolene
metabolites2$X2_alfa.Muurolene <- metabolites2$x1*metabolites$alfa.Muurolene/metabolites$Area.standar.porc
metabolites2$X3_alfa.Muurolene <- metabolites2$X2_alfa.Muurolene*metabolites$microlitres.metabolites/metabolites$inyection.microlitres
metabolites2$X4_alfa.Muurolene <- metabolites2$X3_alfa.Muurolene*100/metabolites$We.needle.grams
metabolites2$X5_alfa.Muurolene <- metabolites2$X4_alfa.Muurolene*1000
    
# Data frame con datos a partir de la [] del estandar 
metabolites3<- data.frame(metabolites2$C.A.T,metabolites2$C.A,metabolites2$Season,metabolites2$Condition,metabolites2$A.expo, metabolites2$Sample, metabolites2$mmet,
                          metabolites2$X4_beta.pinene, metabolites2$X4_L.alfa.bornyl.acetate , metabolites2$X4_beta.Caryophyllene.oxide, metabolites2$X4_alfa.Caryophyllene,
                          metabolites2$X4_beta.Cubebene,metabolites2$X4_alfa.Cubenene,metabolites2$X4_delta.Cadinene,metabolites2$X4_alfa.Muurolene)


# Change col names
colnames(metabolites3)<- c("C-A-T","C-A","Season","Condition", "A.exposition", "Sample","n.GC-MS",
                           "beta.pinene","L.alfa.bornyl.acetate","beta.Caryophyllene.oxide","alfa.Caryophyllene",
                           "beta.Cubebene","alfa.Cubenene","delta.Cadinene","alfa.Muurolene")

# Export data
write.table(metabolites3, "../../metadata/calculate_relative_abs.txt", sep="\t", row.names=FALSE)

