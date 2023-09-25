#!/usr/bin/env bash

for file in *motifsHYP3_v8.fasta;
do
	outfile=$(echo "$file" | awk '{split($0,a,"."); print a[1]"_clipped.fasta"}')
	../motif_clip_HYP3_v8.R $file $outfile
	wait
done



