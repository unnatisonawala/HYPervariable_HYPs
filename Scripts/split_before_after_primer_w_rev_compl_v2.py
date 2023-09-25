#!/usr/bin/env python3 
import csv
import sys 
import pandas as pd
import os 
import re
from Bio.SeqIO.FastaIO import FastaTwoLineParser
from Bio.Seq import Seq
file = sys.argv[1]
life = sys.argv[2]
gene = sys.argv[3]
file = pd.read_csv(file, sep='\t', header = None)
#dir = sys.argv[2]
for i in range(len(file)):
    forward=file.iloc[i,2]
    reverse_rev_comp=file.iloc[i,6]
    reverse=file.iloc[i,4]
    forward_rev_comp=file.iloc[i,5]
    comb1 = [forward, reverse_rev_comp]
    comb2 = [reverse, forward_rev_comp]
    #print(forward + '\n' + reverse_rev_comp + '\n' + reverse + '\n' + forward_rev_comp)
    #print(comb1)
    j= i+1
    k=0
    l=0
    n=0
    o=0
    test_write = open(str(life) + "_" + str(gene) +  "_reads_v3_edited_rev_comp/" + str(life) + "_" + str(gene) + "_edited_%s.fasta"%j, "a") 
    with open(str(life) + "_" + str(gene) + "_reads_v3/" +str(life) + "_" + str(gene) + "_%s.fasta"%j) as handle:
        for title, seq in FastaTwoLineParser(handle):
#             if forward in seq:
#               k= k+1
            if all(x in seq for x in comb1):
                n= n+1
                seq2 = forward + seq.split(forward,1)[1]
                seq2 = seq2.rsplit(reverse_rev_comp, 1)[0] + reverse_rev_comp
            #elif forward_rev_comp in seq:
            #    l=l+1
            elif all(x in seq for x in comb2):
                seq2 = reverse + seq.split(reverse, 1)[1]
                seq2 = seq2.rsplit(forward_rev_comp, 1)[0] + forward_rev_comp
                seq2 = Seq(seq2)
                seq2 = seq2.reverse_complement()
                seq2 = str(seq2)
                o = o+1
             print(title + "\n" + seq2)
             print('k:' + str(k) + 'n:' + str(n) + 'l:' + str(l) + 'o:' + str(o))
             test_write.write(">" + title + "\n" + seq2 + "\n")



