#!/bin/bash

# Example execution: ./correct_fastq.sh
# This script requires RACER, which can be easily installed via conda with:
#   conda install -c mmolnar racer
# Note that a text file called `samples.txt` must exist in the working directory
# and contain the name of all the samples to be corrected, one per line.
# Before using RACER, make sure that the input FastQ files are decompressed.

fastqUnc=fastq_flt
outDir=fastq_corrected

while read sam; do
  racer $fastqUnc/${sam}_R1.fastq outDir/${sam}_R1.fastq 25000
  racer $fastqUnc/${sam}_R2.fastq outDir/${sam}_R2.fastq 25000
  gzip outDir/${sam}_R1.fastq
  gzip outDir/${sam}_R2.fastq
done <samples.txt

