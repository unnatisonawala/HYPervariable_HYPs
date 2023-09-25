#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
library(Biostrings)
library(stringr)
infile <- args[1]
outfile <- args[2]
test <- readDNAStringSet(infile)
print(test[1])
vmatchPattern2 <- function(pattern, subject,
                              max.mismatch=0, min.mismatch=0,
                              with.indels=FALSE, fixed=TRUE,
                              algorithm="auto"){
     if (!is(subject, "XStringSet"))
         subject <- Biostrings:::XStringSet(NULL, subject)
     algo <- Biostrings:::normargAlgorithm(algorithm)
     if (Biostrings:::isCharacterAlgo(algo))
         stop("'subject' must be a single (non-empty) string ",
              "for this algorithm")
     pattern <- Biostrings:::normargPattern(pattern, subject)
     max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
     min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch,
max.mismatch)
     with.indels <- Biostrings:::normargWithIndels(with.indels)
     fixed <- Biostrings:::normargFixed(fixed, subject)
     algo <- Biostrings:::selectAlgo(algo, pattern,
                                     max.mismatch, min.mismatch,
                                     with.indels, fixed)
     C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject,
                     max.mismatch, min.mismatch,
                     with.indels, fixed, algo,
                     "MATCHES_AS_RANGES",
                     PACKAGE="Biostrings")
     unlisted_ans <- IRanges(start=unlist(C_ans[[1L]],
use.names=FALSE),
                             width=unlist(C_ans[[2L]],
use.names=FALSE))
     relist(unlisted_ans, C_ans[[1L]])
   }
#yerggg
motif1 <- vmatchPattern2(pattern = "tatgagcgcggaggcgga", test, min.mismatch = 0, max.mismatch = 2, with.indels = TRUE, fixed = "subject" )
#snrggg or sdrggg 
#edited compared to previous versions to include all possible Gly. Not focussing on nucleotide variants at the moment.
motif2 <- vmatchPattern2(pattern = "agyraycgcggnggnggn", test, max.mismatch = 1, with.indels = TRUE, fixed = "subject")
#sdrgd or rdrgd (b= g/c/t)
motif3<- vmatchPattern2(pattern = "mgygaccgyggagab", test, max.mismatch = 1, with.indels = TRUE,fixed = "subject")
#reggd
motif4<- vmatchPattern2(pattern = "cgtgaaggcggagac", test, max.mismatch = 1, with.indels = TRUE,fixed = "subject")
#rdnkrg
motif5 <- vmatchPattern2(pattern = "cgtgacaataagcgcgga", test, max.mismatch = 2, with.indels = TRUE, fixed = "subject")

#specify motif2 further if the second codon doesn't have errors at 1st pos.
#6 --> SNRGGG
#7 --> SDRGGG
#snrggg
motif6 <-vmatchPattern(pattern = "agyaaycgcggnggnggn", test, max.mismatch = 0, fixed = "subject")
#sdrggg
motif7 <-vmatchPattern(pattern = "agygaycgcggnggnggn", test, max.mismatch = 0, fixed = "subject")


#leaving numbers 8,9,0 to specify motif3
#8--> RDRGD
#9 --> SDRGD
#0 --> SDRGE
#rdrgd
motif8<- vmatchPattern2(pattern = "cgygaccgyggagay", test, max.mismatch = 0, with.indels = TRUE,fixed = "subject")
#sdrgd
motif9<- vmatchPattern2(pattern = "agygaccgyggagay", test, max.mismatch = 0, with.indels = TRUE,fixed = "subject")
#sdrge
motif0<- vmatchPattern2(pattern = "agygaccgyggagar", test, max.mismatch = 0, with.indels = TRUE,fixed = "subject")
motif_list <- list("motif1", "motif3", "motif4", "motif2", "motif5")

for (h in 1:length(test)){
  temp <- as.character(test[h])
  for (j in motif_list){
    motif_number <- str_replace(as.character(j), "motif", "")
    motif_number <- as.integer(motif_number)
    motif_number_2 <- NA
    motif <- get(j)
    for (i in 1:nrow(as.data.frame(motif[h]))){
         if (nrow(as.data.frame(motif[h])) == 0){
           break
         }
         x <- as.data.frame(motif[h])[i,3]
         y <- as.data.frame(motif[h])[i,4]
         z <- as.data.frame(motif[h])[i,5]
         if(motif_number == 2){
           if(str_sub(test[h],(x+3),(x+3)) == "A"){
             motif_number_2 <- 6
           } else if(str_sub(test[h],(x+3),(x+3)) == "G"){
             motif_number_2 <- 7
           }
         }
         if(motif_number == 3){
           
           if(str_sub(test[h],x,x) == "C"){
             motif_number_2 <- 8
           } else if(str_sub(temp,x,x) == "A"){
             
             motif_number_2 <- 9
             if(str_sub(test[h],(x+14),(x+14)) == "A" || str_sub(test[h],(x+14),(x+14)) == "G"  ){
               motif_number_2 <- 0 
             }
           }
           
         }
         
         if(!is.na(motif_number_2)){
           str_sub(temp, x, y) <- strrep(as.character(motif_number_2), z)
           motif_number_2 <- NA
         } else {
         str_sub(temp, x, y) <- strrep(as.character(motif_number), z) }
        }
  }
  cat(">",names(test[h]),"\n", sep='', file = outfile, append = TRUE)
  cat(temp, "\n", sep='',file= outfile, append = TRUE)
}
