# Example execution: This script must be run inside an R session to explore the
# quality of the MiSeq (deep amplicon sequencing) paired-end reads. Just change 
# the value of variables `fastqDir` and `outputDir`. We recommend running this
# script on ~12 samples at a time (i.e. put the forward and reverse FastQ files
# of ~12 samples in the `fastqDir` folder).

library(dada2)
library(ggplot2)

fastqDir <- "fastq/"
outputDir <- "fastq_flt/"

# Forward and reverse fastq filenames
fnFs <- sort(list.files(fastqDir, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(fastqDir, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(fnFs) + geom_vline(xintercept = 150)
plotQualityProfile(fnRs) + geom_vline(xintercept = 150)

# Filter and trim

# Place filtered files in output directory
filtFs <- file.path(outputDir, paste0(sample.names, "_flt_R1.fastq.gz"))
filtRs <- file.path(outputDir, paste0(sample.names, "_flt_R2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Truncate at 150 (where quality drops in both F and R), trim the first 5 bases, 
# filter out Ns, accept a maximum of 5 expected errors, truncate reads at 
# the first instance of a quality score less than 2
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     trimLeft = 5, maxN=0, maxEE=c(5,5), truncQ=2, 
                     rm.phix=FALSE, compress=TRUE, multithread=TRUE)
out

# Now check the quality of the filtered reads
fnFs_flt <- sort(list.files(outputDir, pattern="_R1.fastq", full.names = TRUE))
fnRs_flt <- sort(list.files(outputDir, pattern="_R2.fastq", full.names = TRUE))
plotQualityProfile(fnFs_flt) + geom_vline(xintercept = 145)
plotQualityProfile(fnRs_flt) + geom_vline(xintercept = 145)
