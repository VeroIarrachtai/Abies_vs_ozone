### This script parse and process html files to dataframe ###
## September 2019
## Verónica Reyes Galindo/Cristian 

## Load libraries
library(tidyr)
library(dplyr)
library(stringr)

## Defines the name of the function, which accepts as parameters: 
## file = "../file.txt", cas = is the number of repeats that you want (example, cas = 3) 

parse_htm_file <- function(file, cas){
  
  rawHTML <- readLines(file)   ## Read text lines from file
  
  rt_grep <- grep("[0-9]+[.][0-9]+[ ][C]", rawHTML) ## Search for matches to pattern to extract Pk#, RT and Area% 
  
  rt_ext <- rawHTML[rt_grep] ## create a list with the matches 
  
  list_rt <- strsplit(rt_ext, "\\s+") ## split each object within the list 
  
  df_rt <- data.frame() ## create a dataframe empty to use in the loop
  
  for (i in 1:length(list_rt)){  ## get the length of the list
    
    Pk <- as.numeric(list_rt[[i]][2]) ## get the value for Pk#
    
    RT <- as.numeric(list_rt[[i]][3])  ## get the value for RT
    
    Area <- as.numeric(list_rt[[i]][4]) ## get the value for Area%
    
    df_first <- data.frame(cbind(Pk, RT, Area)) ##  combine Pk#, RT and Area% by columns 
    
    df_rt <- rbind(df_rt, df_first) ##  combine by row the df_first in df_rt
    
  }
  
  df_rt <- df_rt[rep(seq_len(nrow(df_rt)), each=cas),] ## RepeatPk#, RT and Area% the number of times given by cas
  
  rownames(df_rt) <- NULL ## delete row names
  
  cas_grep <- grep("[0-9]+[ ]+[0-9]+[-][0-9]+[-][0-9]+[ ]+[0-9]+", rawHTML, value = TRUE)
  
  list_cas <- str_extract(cas_grep, "[0-9]+[ ]+[0-9]+[-][0-9]+[-][0-9]+[ ]+[0-9]+")
  
  new <- gsub("\\s+", " ", list_cas) ## Delete more than one space
  
  new <- trimws(new) ## Remove leading and/or trailing whitespace from character strings.
  
  new <- gsub(" +", "_", new) ## replace space with _
  
  tab_cas <- as.data.frame(new) 
  
  tab_cas <- separate(tab_cas, new, c("Ref", "CAS", "Qual"), "_") ## separate Ref#, CAS# and Qual into columns 
  
  tab_cas$Ref <- as.numeric(tab_cas$Ref) 
  
  tab_cas$Qual <- as.numeric(tab_cas$Qual)
  
  final_df <- cbind(df_rt, tab_cas) ## combine df_rt (Pk#, RT and Area%) and  tab_cas (Ref#, CAS# and Qual) by columns
  
  return(final_df) 
  
}
