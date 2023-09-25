#!/usr/bin/env bash

for file in *;
do
	outfile=$(echo "$file" | awk '{split($0,a,"."); print a[1]"_motifsHYP3_v8.fasta"}')
	../motif_analysis_HYP3_v8.R $file $outfile
	wait
done



