# Example execution: Rscript trim_and_filter.R --dir fastq --out fastq_flt
# Note that this script trims and filters all the FastQ files in the given 
# directory and also outputs sample statistics in the working directory.

library(dada2)
library(optparse)

option_list <- list(
  make_option(
    c("-d", "--dir"),
    type = "character",
    default = ".",
    help = "The directory with the FastQ files that will be trimmed and filtered [default = current directory]."
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = "fastq_flt",
    help = "The directory where trimmed FastQ files will be written [default = fastq_flt]."
  )
)

opt_parser <-
  OptionParser(usage = "Rscript trim_and_filter.R [options]",
               option_list = option_list)
opt <-  parse_args(opt_parser)

fastq_path <- opt$dir
out_path <- opt$out

# Forward and reverse fastq filenames
fnFs <-
  sort(list.files(fastq_path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <-
  sort(list.files(fastq_path, pattern = "_R2.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# Form the names of the output files
if (!dir.exists(out_path)) {
  dir.create(out_path)
}
filtFs <- file.path(out_path, paste0(sample.names, "_R1.fastq"))
filtRs <- file.path(out_path, paste0(sample.names, "_R2.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Trim and filter
out <-
  filterAndTrim(
    fnFs,
    filtFs,
    fnRs,
    filtRs,
    truncLen = c(150, 150),
    trimLeft = 5,
    maxN = 0,
    maxEE = c(5, 5),
    truncQ = 2,
    rm.phix = FALSE,
    compress = FALSE,
    multithread = TRUE
  )

# Save cleaning stats in the current directory
out_stats <- tibble::tibble(sample_name = sample.names,
                    reads_in = out[, "reads.in"],
                    reads_out = out[, "reads.out"])
readr::write_csv(out_stats, "./sample_stats.csv")
