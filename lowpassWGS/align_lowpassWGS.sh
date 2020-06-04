#!/bin/bash

# Example execution: ./align_lowpassWGS.sh L_C12.01

# Variable definitions
sample=$1

fastqDir=/home/user_name/lowpassWGS

# Genome downloaded from the Illumina iGenomes repository
genome=genome
genomeDir=Hsapiens_UCSC_hg19/Sequence/BWAIndex

# Output directory
resDir=/home/user_name/lowpassWGS/bams

bwa mem $genomeDir/${genome}.fa $fastqDir/${sample}.fastq.gz > $resDir/${sample}.sam

# Convert to BAM
samtools view -S -b $resDir/${sample}.sam > $resDir/${sample}.bam

# Sort
samtools sort $resDir/${sample}.bam -o $resDir/${sample}_sorted.bam

# Generate an index
samtools index $resDir/${sample}_sorted.bam
