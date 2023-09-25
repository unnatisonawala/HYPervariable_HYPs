#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
print(args)
infile <- args[1]
outfile <- args[2]
library(Biostrings)
library(stringr)
# You need stringr version 1.4.1 so that it also takes in DNA stringsets as characters. The newer version 1.5 throws an error. Download older version using package remotes and then install_version
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

test <- readDNAStringSet(infile)
#test <- readDNAStringSet("/Users/endurance_v3/OneDrive - University of Cambridge/Data/G.pall/J2_HYP3_edited_53.fasta")
#SDEKEKKYG
#1. for S 
motif1 <- vmatchPattern2(pattern = "AGCGATGArAAGGAAAAAAAGTACGGT", test, min.mismatch = 0, max.mismatch = 1, with.indels = TRUE, fixed = "subject" )

#_DEK__EKKYG
#2. for N
#3. for K
#4. for S
motif2 <- vmatchPattern2(pattern = "armgatvrnaarrarrargaraaraartacggt", test, min.mismatch = 0, max.mismatch = 1, with.indels = TRUE, fixed = "subject" )
#motif2 <- vmatchPattern2(pattern = "armgatgrdaaggaarargagaagaartacggt", test, min.mismatch = 0, max.mismatch = 1, with.indels = TRUE, fixed = "subject" )
#5 --> SDD
#6 --> SDQ
#motif5 <- vmatchPattern2(pattern = "agcgatsakaargaagaggagaagaaatacggt", test, min.mismatch = 0, max.mismatch = 0, with.indels = TRUE, fixed = "subject" )

#NEKEEEKKYG
motif7 <- vmatchPattern2(pattern = "aatgagaaggaagaggagaagaaatacggt", test, min.mismatch = 0, max.mismatch = 1, with.indels = TRUE, fixed = "subject" )

#NDGKEKEKKYG
#motif8 <- vmatchPattern2(pattern = "AACGATGGNAAGGAAAAGGAAAAGAAGTACGGT", test, min.mismatch = 0, max.mismatch = 0, with.indels = TRUE, fixed = "subject" )

#SDENEJEESGK....
motif9 <- vmatchPattern2(pattern = "agcgatgagaatgaaaaggaggaaagtggtaaaggatccggaggt", test, min.mismatch = 0, max.mismatch = 1, with.indels = TRUE, fixed = "subject" )

motif_list <- list("motif1", "motif7", "motif2", "motif9")

motif_number_2 <- NA
for (h in 1:length(test)){
  temp <- as.character(test[h])
  for (j in motif_list){
    motif_number <- str_replace(as.character(j), "motif", "")
    motif_number <- as.integer(motif_number)
    motif <- get(j)
    for (i in 1:nrow(as.data.frame(motif[h]))){
         if (nrow(as.data.frame(motif[h])) == 0){
           break
         }
         x <- as.data.frame(motif[h])[i,3]
         y <- as.data.frame(motif[h])[i,4]
         z <- as.data.frame(motif[h])[i,5]
         
         if(motif_number == 2){
    if(str_sub(test[h],(x+2),(x+2)) == "A" || str_sub(test[h],(x+2),(x+2)) == "G"){
             motif_number_2 <- 3 #KDEK
           }
    else if (str_sub(test[h],(x+1),(x+1)) == "G"){
             if (str_sub(test[h],(x+6),(x+6)) == "C" ){
                motif_number_2 <- 6 #SDQ
           }
             else if (str_sub(test[h],(x+6),(x+6)) == "G"){
                if (str_sub(test[h],(x+8),(x+8)) == "C" || str_sub(test[h],(x+8),(x+8)) == "T"){
                    motif_number_2 <- 5 #SDD
                }
                else if (str_sub(test[h],(x+8),(x+8)) == "A" || str_sub(test[h],(x+8),(x+8)) == "G"){
                    motif_number_2 <- 4 #SDE
                }

             }
              else if (str_sub(test[h],(x+6),(x+6)) == "A" ){
                motif_number_2 <- 0
              }
      
           }
    else if(str_sub(test[h],(x+2),(x+2)) == "C" || str_sub(test[h],(x+2),(x+2)) == "T"){
             if (str_sub(test[h],(x+7),(x+7)) == "G" ){
                motif_number_2 <- 8
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
