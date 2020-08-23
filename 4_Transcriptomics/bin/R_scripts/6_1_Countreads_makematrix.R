## Data to diferantial expresion (mapping with BWA)
## 15 Oct 2018
## Veronica Reyes


#Load libraries

library(dplyr)
library(readr)

# Load files

TC01_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR1_S179_sw10L50.genesorder.txt", header= FALSE)
TC02_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR2_S180_sw10L50.genesorder.txt", header= FALSE)
TC03_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR3_S181_sw10L50.genesorder.txt", header= FALSE)
TC04_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR4_S182_sw10L50.genesorder.txt", header= FALSE)
TC05_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR5_S183_sw10L50.genesorder.txt", header= FALSE)

TS01_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR11_S189_sw10L50.genesorder.txt", header= FALSE)
TS02_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR12_S190_sw10L50.genesorder.txt", header= FALSE)
TS05_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR13_S191_sw10L50.genesorder.txt", header= FALSE)

DC01_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR6_S184_sw10L50.genesorder.txt", header= FALSE)
DC02_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR7_S185_sw10L50.genesorder.txt", header= FALSE)
DC03_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR8_S186_sw10L50.genesorder.txt", header= FALSE)
DC04_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR9_S187_sw10L50.genesorder.txt", header= FALSE)
DC05_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR10_S188_sw10L50.genesorder.txt", header= FALSE)

DS01_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR14_S192_sw10L50.genesorder.txt", header= FALSE)
DS02_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR15_S193_sw10L50.genesorder.txt", header= FALSE)
DS04_15<-read.delim("../../data/TXT/GENES_ORDER/DPVR16_S194_sw10L50.genesorder.txt", header= FALSE)

TC01_17<-read.delim("../../data/TXT/GENES_ORDER/DPVR17_S195_sw10L50.genesorder.txt", header= FALSE)
DC04_17<-read.delim("../../data/TXT/GENES_ORDER/DPVR18_S196_sw10L50.genesorder.txt", header= FALSE)

head(DC04_17)

# Change columns name´s
colnames(DC01_15)<- c("number_reads", "name_gen")
colnames(DC02_15)<- c("number_reads", "name_gen")
colnames(DC03_15)<- c("number_reads", "name_gen")
colnames(DC04_15)<- c("number_reads", "name_gen")
colnames(DC05_15)<- c("number_reads", "name_gen")

colnames(DS01_15)<- c("number_reads", "name_gen")
colnames(DS02_15)<- c("number_reads", "name_gen")
colnames(DS04_15)<- c("number_reads", "name_gen")

colnames(TC01_15)<- c("number_reads", "name_gen")
colnames(TC02_15)<- c("number_reads", "name_gen")
colnames(TC03_15)<- c("number_reads", "name_gen")
colnames(TC04_15)<- c("number_reads", "name_gen")
colnames(TC05_15)<- c("number_reads", "name_gen")

colnames(TS01_15)<- c("number_reads", "name_gen")
colnames(TS02_15)<- c("number_reads", "name_gen")
colnames(TS05_15)<- c("number_reads", "name_gen")

colnames(TC01_17)<- c("number_reads", "name_gen")
colnames(DC04_17)<- c("number_reads", "name_gen")
head(DC04_17)

# Changes columns order

all<- data.frame(DC01_15$name_gen, DC01_15$number_reads)
head(all)
all$DC_2 <- DC02_15$number_reads[match(all$DC01_15.name_gen, DC02_15$name_gen)]
all$DC_3 <- DC03_15$number_reads[match(all$DC01_15.name_gen, DC03_15$name_gen)]
all$DC_4 <- DC04_15$number_reads[match(all$DC01_15.name_gen, DC04_15$name_gen)]
all$DC_5 <- DC05_15$number_reads[match(all$DC01_15.name_gen, DC05_15$name_gen)]

all$DS_1 <- DS01_15$number_reads[match(all$DC01_15.name_gen, DS01_15$name_gen)]
all$DS_2 <- DS02_15$number_reads[match(all$DC01_15.name_gen, DS02_15$name_gen)]
all$DS_4 <- DS04_15$number_reads[match(all$DC01_15.name_gen, DS04_15$name_gen)]

all$TC_1 <- TC01_15$number_reads[match(all$DC01_15.name_gen, TC01_15$name_gen)]
all$TC_2 <- TC02_15$number_reads[match(all$DC01_15.name_gen, TC02_15$name_gen)]
all$TC_3 <- TC03_15$number_reads[match(all$DC01_15.name_gen, TC03_15$name_gen)]
all$TC_4 <- TC04_15$number_reads[match(all$DC01_15.name_gen, TC04_15$name_gen)]
all$TC_5 <- TC05_15$number_reads[match(all$DC01_15.name_gen, TC05_15$name_gen)]

all$TS_1 <- TS01_15$number_reads[match(all$DC01_15.name_gen, TS01_15$name_gen)]
all$TS_2 <- TS02_15$number_reads[match(all$DC01_15.name_gen, TS02_15$name_gen)]
all$TS_5 <- TS05_15$number_reads[match(all$DC01_15.name_gen, TS05_15$name_gen)]

all$TC17 <- TC01_17$number_reads[match(all$DC01_15.name_gen, TC01_17$name_gen)]
all$DC47 <- DC04_17$number_reads[match(all$DC01_15.name_gen, DC04_17$name_gen)]

head(all)
colnames(all)[1]<- ""
colnames(all)[2]<- "DC_1"
head(all)

# Change NA to 0´s

all<- as.data.frame(all, na.rm=TRUE)
all[is.na(all)]<- 0

#Remove no match
all<- all[-(1),]

# Count reads per sample

sum(all$DC_1)
sum(all$DC_2)
sum(all$DC_3)
sum(all$DC_4)
sum(all$DC_5)

sum(all$DS_1)
sum(all$DS_2)
sum(all$DS_4)

sum(all$TC_1)
sum(all$TC_2)
sum(all$TC_3)
sum(all$TC_4)
sum(all$TC_5)

sum(all$TS_1)
sum(all$TS_2)
sum(all$TS_5)

# Count genes with 0 reads
nrow(all[all$DC_1 == 0,])
nrow(all[all$DC_2 == 0,])
nrow(all[all$DC_3 == 0,])
nrow(all[all$DC_4 == 0,])
nrow(all[all$DC_5 == 0,])

nrow(all[all$DS_1 == 0,])
nrow(all[all$DS_2 == 0,])
nrow(all[all$DS_4 == 0,])

nrow(all[all$TC_1 == 0,])
nrow(all[all$TC_2 == 0,])
nrow(all[all$TC_3 == 0,])
nrow(all[all$TC_4 == 0,])
nrow(all[all$TC_5 == 0,])

nrow(all[all$TS_1 == 0,])
nrow(all[all$TS_2 == 0,])
nrow(all[all$TS_5 == 0,])

# Export table to .txt format
write.table(all, file="../../data/allreadsgenes.txt", sep = "\t", eol = "\n", dec = ".",
            row.names = TRUE, col.names = TRUE)
