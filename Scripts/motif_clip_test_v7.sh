#!/usr/bin/env bash

for file in *_motifs_v3.fasta;
do
	outfile=$(echo "$file" | awk '{split($0,a,"."); print a[1]"v7_clipped.fasta"}')
	../motif_clip_v7.R $file $outfile
	wait
done



