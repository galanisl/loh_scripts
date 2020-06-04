# Example execution: This script must be run inside an R session to generate
# a copy-number profile for each sample. Just change the value of variable 
# `sample_name` and the location of the corresponding BAM file if necessary.

library(QDNAseq)

# Name of the sample to be processed
sample_name <- "L_C12.01"
# Location of the sorted BAM file for the sample
bam <- paste0("bams/", sample_name, "_sorted.bam")

bins <- getBinAnnotations(binSize = 100)

readCounts <- binReadCounts(bins, bamfiles = bam)
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, 
                                   blacklist = TRUE, chromosomes = c("MT"))
readCountsFiltered <- estimateCorrection(readCountsFiltered)

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

# Plot the full copy number profile
plot(copyNumbersSegmented)

# Focus on chromosomes 5, 6 and 7
chr5 <- chromosomes(copyNumbersSegmented) == 5
chr6 <- chromosomes(copyNumbersSegmented) == 6
chr7 <- chromosomes(copyNumbersSegmented) == 7

plot(copyNumbersSegmented[chr5, ])
plot(copyNumbersSegmented[chr6, ])
plot(copyNumbersSegmented[chr7, ])