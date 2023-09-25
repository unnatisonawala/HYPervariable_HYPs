#!/usr/bin/env bash

for file in *;
do
	outfile=$(echo "$file" | awk '{split($0,a,"."); print a[1]"_motifs_v3.fasta"}')
	../motif_analysis_v3.R $file $outfile
	wait
done



