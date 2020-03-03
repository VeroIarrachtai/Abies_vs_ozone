# Barplot abundance per sample
# Septiembre 2019
# Verónica Reyes

# Load libraries
library(ggplot2)
library(reshape2)

# Load data

metabolites<-read.delim("../../metadata/calculate_relative_abs.txt")

# Plotear


metabolites_mean <- aggregate(metabolites[,8:15], by=list(Factors=metabolites$C.A.T), FUN=mean)
metabolites_sd <- aggregate(metabolites[,8:15], by=list(Factors=metabolites$C.A.T), FUN=sd)
df_mean <- melt(metabolites_mean, id.vars=c("Factors"), variable.name = "metabolite", value.name="value")
df_sd <- melt(metabolites_sd, id.vars=c("Factors"), variable.name = "metabolite", value.name="value")


#SS
levels(df_mean$Factors)
levels(df_mean$Factors)<- c("Contingency Damaged 2015","Moderate concentration Damaged 2015",
                             "Contingency Damaged 2016","Moderate concentration Damaged 2016",
                             "Contingency Tolerant 2015","Moderate concentration Tolerant 2015",
                             "Contingency Tolerant 2016","Moderate concentration Tolerant 2016")


metabolites4<-data.frame(metabolites)[1:20,]

metabolites_mean <- aggregate(metabolites4[,8:15], by=list(Factors=metabolites4$C.A.T), FUN=mean)
metabolites_sd <- aggregate(metabolites4[,8:15], by=list(Factors=metabolites4$C.A.T), FUN=sd)
df_mean <- melt(metabolites_mean, id.vars=c("Factors"), variable.name = "metabolite", value.name="value")
df_sd <- melt(metabolites_sd, id.vars=c("Factors"), variable.name = "metabolite", value.name="value")


limits <- aes(ymax = df_mean[,"value"] + df_sd[,"value"], ymin=df_mean[,"value"] - df_sd[,"value"])


df_mean$labelss <- c("M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016",
                    "M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016",
                    "M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016",
                    "M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016",
                    "M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016",
                    "M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016",
                    "M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016",
                    "M. concentration Damaged 2015","M. concentration Damaged 2016",
                    "M. concentration Tolerants 2015","M. concentration Tolerants 2016")
dev.off()
p <- ggplot(df_mean, aes(metabolite, value, fill = Factors)) +
  geom_bar(position="dodge", stat="identity") + geom_errorbar(limits, position="dodge")
print(p)
p + coord_flip() + facet_wrap(~ Factors) +
  scale_fill_manual(values= c( "#c6003a", "#e98382", "#00901e","#b1e787"))+
  facet_wrap(~ labelss)+
  scale_x_discrete (labels = c('beta.pinene' = expression(beta~'-Pinene'),
                               'L.alfa.bornyl.acetate' = expression('L-'~ alpha ~'-Bornyl acetate'),
                               'beta.Caryophyllene.oxide'= expression(beta~'-Caryophyllene oxide'),
                               'alfa.Caryophyllene' = expression(alpha~'-Caryophyllene'),
                               'beta.Cubebene'= expression(beta~'-Cubebene'),
                               'alfa.Cubenene'= expression(alpha~'-Cubenene'),
                               'delta.Cadinene' = expression(delta~'-Cadinene'),
                               'alfa.Muurolene' = expression(alpha~'-Muurolene')))+
  labs(x="metabolite",y="g/100g Tissue")

ggsave("../../outputs/4.1_barplot_images_SS.png")

#Contingency

metabolites4<-data.frame(metabolites)[21:40,]

metabolites_mean <- aggregate(metabolites4[,8:15], by=list(Factors=metabolites4$C.A.T), FUN=mean)
metabolites_sd <- aggregate(metabolites4[,8:15], by=list(Factors=metabolites4$C.A.T), FUN=sd)
df_mean <- melt(metabolites_mean, id.vars=c("Factors"), variable.name = "metabolite", value.name="value")
df_sd <- melt(metabolites_sd, id.vars=c("Factors"), variable.name = "metabolite", value.name="value")

df_mean$labelss <- c("Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016",
                     "Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016",
                     "Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016",
                     "Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016",
                     "Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016",
                     "Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016",
                     "Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016",
                     "Contingency Damaged 2015","Contingency Damaged 2016",
                     "Contingency Tolerants 2015","Contingency Tolerants 2016")

p <- ggplot(df_mean, aes(metabolite, value, fill = Factors)) +
  geom_bar(position="dodge", stat="identity") + geom_errorbar(limits, position="dodge")
print(p)
p + coord_flip() + facet_wrap(~ Factors) +
  scale_fill_manual(values= c( "#c6003a", "#e98382", "#00901e","#b1e787"))+
  facet_wrap(~ labelss)+
  scale_x_discrete (labels = c('beta.pinene' = expression(beta~'-Pinene'),
                               'L.alfa.bornyl.acetate' = expression('L-'~ alpha ~'-Bornyl acetate'),
                               'beta.Caryophyllene.oxide'= expression(beta~'-Caryophyllene oxide'),
                               'alfa.Caryophyllene' = expression(alpha~'-Caryophyllene'),
                               'beta.Cubebene'= expression(beta~'-Cubebene'),
                               'alfa.Cubenene'= expression(alpha~'-Cubenene'),
                               'delta.Cadinene' = expression(delta~'-Cadinene'),
                               'alfa.Muurolene' = expression(alpha~'-Muurolene')))+
  labs(x="metabolite",y="g/100g Tissue")

ggsave("../../outputs/4.1_barplot_images_conti.png")

