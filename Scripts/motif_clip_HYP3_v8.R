#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
library(Biostrings)
library(stringr)
motif_list <- list("motif1", "motif3", "motif4", "motif2", "motif5","motif6","motif7","motif8","motif9","motif0")
test<- readBStringSet(infile)
for (h in 1:length(test)){
  temp <- as.character(test[h])
  for (x in 0:(length(motif_list)-1)){
        for (i in 132:22){
                if (x == 2 || x == 3 || x == 4 || x == 5 || x == 6 || x == 8 || x==0){
                        temp <- str_replace_all(temp, str_dup(x,i), str_dup(x, round(i/33)) )
                }
                else if (x == 1){
                        temp <- str_replace_all(temp, str_dup(x,i), str_dup(x, round(i/27)) )
                }
                else if (x == 7){
                        temp <- str_replace_all(temp, str_dup(x,i), str_dup(x, round(i/30)) )
                }
                else if (x == 9){

                         temp <- str_replace_all(temp, str_dup(x,i), str_dup(x, round(i/45)) )
                }
        }
  }
  temp <- gsub("(?<=\\d)[[:alpha:]]{24,36}(?=\\d)","_",temp, perl=TRUE)
  temp <- gsub("(?<=\\d)[[:alpha:]]{50,70}(?=\\d)","__",temp, perl=TRUE)
  temp <- gsub("(?<=\\d)[[:alpha:]]{90,130}(?=\\d)","___",temp, perl=TRUE)
  temp <- gsub("[[:alpha:]]","",temp)
  cat(temp,"\n",sep="", file= outfile ,append=TRUE)
}
