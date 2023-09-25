## Splitting and cleaning up barcoded Pacbio reads into fasta files for individuals
* The entire sequencing dataset was mapped to Newton assembly 
```bash
minimap2 -ax map-hifi /path/to/file/Gpal_newton_newton.scaffolds.fa m64101_220402_211610.hifi_reads.fasta
```
* The bam file was sorted and indexed 
* Reads mapping to HYP1 and HYP3 loci were filtered using samtools view function to filter out any reads from off-target amplification. 
* The reads were split by barcodes using the following 
<mark> folder with text files for barcoded primers </mark>
```bash 
./barcode_splittter_v3.sh {text file with primers specific to HYP1 or HYP3 and life stage} {data/fasta file with reads} {directory to save to} {life stage} {gene - HYP1 or HYP3}
```
* The following python script was used to clip sequences before the forward or after the reverse primer. Reverse strand reads were reverse complemented.
```bash 
./split_before_after_primer_w_rev_compl_v2.py {text file with primers specific to HYP1 or HYP3 and life stage} {life stage} {gene -HYP1 or HYP3}

```
'
## Notes on getting HVD as a form of a string of numbers from Pacbio or Nanopore reads
 **For getting HVD for Pacbio reads from individual nematodes**
* Run **motif_analysis_test_v3.sh**, which will then run the Rscript **motif_analsysis_v3.R** on each fasta containing reads from individual nematodes. 
* This will return a new fasta for each where the HVD has now been replaced by numbers. Each number represents a motif and the same number is repeated for the length of the motif to retain the original fasta length. 
* On these "motif annotated" fasta run **motif_clip_test_v7.sh** which in turn will run an Rscript **motif_clip_v7.R** on each file. This clips the fasta to remove the non-HVD portions and also reduce the string of number for each motif down to one. 
***Final file:*** The result should be a file with a string of numbers representing HVD for each read per line. 
N.B: Run motif_analysis_HYP3_v4.sh <mark> and the clipping file </mark> for HYP3
