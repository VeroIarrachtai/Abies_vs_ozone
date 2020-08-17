### This script make a dataframe with all samples of GC-SM analisys###
## September 2019
## Verónica Reyes Galindo/Cristian

## Call the function 
source("Load_html_files.R")

## List the directories 
folders <- list.files("../../data/GC-MS/")

## Create a dataframe empty
htm_df <- data.frame()

for (i in folders){
  ## Since the directories end in _outfiles we have to stay alone with the folder identifier in this case clust_ [0-9]+_ [0-9]+
  name <- stringr::str_extract(string = i, pattern = "LibSrch_[0-9]+")
  
  ## Extract the desired tree, for this we use paste0 to enter the folders where are the trees and the mean as dataframe
  tab_htm <- parse_htm_file(paste0("../../data/GC-MS/", i), 1 )
  
  ## Create a dataframe repeating the identifier according to with the number of row of the table 
  tab_names <- data.frame(rep(name, times = nrow(tab_htm)))
  
  ## Change the name of identifier column
  colnames(tab_names) <- "Id"
  
  ## combine the data frames by columns
  c_htm <- cbind(tab_names, tab_htm)
  
  ## combine the data frames by rows    
  htm_df <- rbind(htm_df, c_htm )
}

write.table(htm_df, "../../metadata/htm_df_allsamples.txt", sep="\t", row.names=FALSE)




