#!/usr/bin/env bash

#run script witthin the directory containing your file with barcode names and fasta file
#use absolute paths for save directory if not under your directory you run script from

file_to_parse=$1
data=$2
directory_to_save=$3
sample=$4
gene=$5
#worm_number=$(cat $file_to_parse | cut -f1 -)
#forward_primer=$(cat $file_to_parse | cut -f3 -)
#reverse_primer=$(cat $file_to_parse | cut -f5 -) 
i=1
while IFS=$'\t' read -r -a line;
do 
	echo "${line[2]}"
	echo "${line[6]}"
	echo "${line[4]}"
	echo "${line[5]}"
	grep -B 1 --no-group-separator "${line[2]}" $data | grep -B 1 --no-group-separator "${line[6]}" - > $4_temp_$i.txt 
	
	grep -B 1 --no-group-separator "${line[4]}" $data | grep -B 1 --no-group-separator "${line[5]}" - > $4_temp_rev_$i.txt
	cat $4_temp*$i.txt > $3/$4_$5_$i.fasta
	rm $4_temp*$i.txt
  	#grep -B 1 --no-group-separator "${line[6]}$" $data > test_temp_$i.txt 
	i=$((i+1))
done < "$file_to_parse"



