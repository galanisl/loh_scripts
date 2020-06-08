#!/bin/bash

# Example execution: ./call_snps.sh G_C16.05

# Variable definitions
sample=$1

fastqDir=/home/user_name/working_dir/fastq_corrected

# Genome downloaded from the Illumina iGenomes repository
genome=genome
idxDir=/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex
genomeDir=/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta

bamDir=/home/user_name/working_dir/bams
snpDir=/home/user_name/working_dir/vcf
outDir=/home/user_name/working_dir/vcf_filtered

# Alignment
bwa mem -t 4 $idxDir/${genome}.fa \
  $fastqDir/${sample}_R1.fastq.gz $fastqDir/${sample}_R2.fastq.gz > tmp/${sample}.sam

# Convert to BAM, sort by position and index
samtools sort -O bam tmp/${sample}.sam -o tmp/${sample}_firstsort.bam
samtools index tmp/${sample}_firstsort.bam 

# Focus on region of interest
samtools view -h tmp/${sample}_firstsort.bam \
  chr3:121584000-121586000 chr6:31156000-31181000 -o tmp/${sample}_zoom.bam

# Sort by name, clean up read pairing information and flags, sort by position
samtools sort -n tmp/${sample}_zoom.bam | \
  samtools fixmate -m - - | \
  samtools sort - -o $bamDir/${sample}.bam

# Generate an index
samtools index $bamDir/${sample}.bam

# FastQC the bam file
fastqc $bamDir/${sample}.bam --outdir rep_fastqc

# SNP calling
bcftools mpileup --threads 16 --max-depth 2000 -a'AD,DP,ADF,ADR' \
  -Ou -f $genomeDir/${genome}.fa $bamDir/${sample}.bam | \
  bcftools call --threads 16 -mv -V 'indels' -Ov -o $snpDir/${sample}.vcf

# Filter based on reads necessary for 5x coverage and a MQ of 50
vcfutils.pl varFilter -Q 50 -d 10 $snpDir/${sample}.vcf > $outDir/${sample}_filtered.vcf

