# Example execution: This script must be run inside an R session to generate
# a digital karyotype for all samples in the gene expression matrix (CSV file
# read in line 10)

library(dplyr)
library(readr)
library(annotables)
library(pheatmap)

rna <- read_csv("GSE100118_scRNA_pou5f1crispr_rpkm_170603.csv")

rna_expr <- rna %>% 
  select(`2.1`:`C24.5`)

rna_genes <- rna %>% 
  select(Ensembl, Gene)

# Go from FPKM to TPM
fpkm_sum <- colSums(rna_expr)
tpm <- (as.matrix(rna_expr) %*% diag(1/fpkm_sum)) * 10^6
colnames(tpm) <- colnames(rna_expr)
rownames(tpm) <- rna_genes$Ensembl

# Use annotables data to get chromosome location of each gene
grch38_flt <- grch38 %>% 
  select(ensgene, chr, start, end)
rna_genes <- left_join(rna_genes, grch38_flt, by = c("Ensembl" = "ensgene"))
rna_genes <- rna_genes %>% 
  filter(!duplicated(rna_genes))

# Remove no-show genes, genes without chromosome info, mitochondrial and PAR genes
no_show <- rowSums(tpm) == 0
no_chr <- is.na(rna_genes$chr)
mito <- rna_genes$chr == "MT"
par_genes <- rna_genes$Gene %in% c("PLCXD1","GTPBP6","PPP2R3B","SHOX","CRLF2",
                                   "CSF2RA","IL3RA","SLC25A6","ASMTL","P2RY8",
                                   "CXYorf3","ASMT","DHRSXY","ZBED1","CD99","XG")

tpm <- tpm[!(no_show | no_chr | mito | par_genes), ]
rna_genes <- rna_genes[!(no_show | no_chr | mito | par_genes), ]

# Add chromosome arm information to each gene
rna_genes <- rna_genes %>% 
  mutate(arm = case_when(
    chr == "1" & end < 133000000 ~ "p", chr == "1" & end > 133000000 ~ "q",
    chr == "2" & end < 92200000 ~ "p", chr == "2" & end > 92200000 ~ "q",
    chr == "3" & end < 91000000 ~ "p", chr == "3" & end > 91000000 ~ "q",
    chr == "4" & end < 50500000 ~ "p", chr == "4" & end > 50500000 ~ "q",
    chr == "5" & end < 48000000 ~ "p", chr == "5" & end > 48000000 ~ "q",
    chr == "6" & end < 59800000 ~ "p", chr == "6" & end > 59800000 ~ "q",
    chr == "7" & end < 60000000 ~ "p", chr == "7" & end > 60000000 ~ "q",
    chr == "8" & end < 45000000 ~ "p", chr == "8" & end > 45000000 ~ "q",
    chr == "9" & end < 48000000 ~ "p", chr == "9" & end > 48000000 ~ "q",
    chr == "10" & end < 40000000 ~ "p", chr == "10" & end > 40000000 ~ "q",
    chr == "11" & end < 52500000 ~ "p", chr == "11" & end > 52500000 ~ "q",
    chr == "12" & end < 36200000 ~ "p", chr == "12" & end > 36200000 ~ "q",
    chr == "13" & end < 16000000 ~ "p", chr == "13" & end > 16000000 ~ "q",
    chr == "14" & end < 16000000 ~ "p", chr == "14" & end > 16000000 ~ "q",
    chr == "15" & end < 18000000 ~ "p", chr == "15" & end > 18000000 ~ "q",
    chr == "16" & end < 40000000 ~ "p", chr == "16" & end > 40000000 ~ "q",
    chr == "17" & end < 25000000 ~ "p", chr == "17" & end > 25000000 ~ "q",
    chr == "18" & end < 18000000 ~ "p", chr == "18" & end > 18000000 ~ "q",
    chr == "19" & end < 25000000 ~ "p", chr == "19" & end > 25000000 ~ "q",
    chr == "20" & end < 27000000 ~ "p", chr == "20" & end > 27000000 ~ "q",
    chr == "21" & end < 11000000 ~ "p", chr == "21" & end > 11000000 ~ "q",
    chr == "22" & end < 14000000 ~ "p", chr == "22" & end > 14000000 ~ "q",
    chr == "X" & end < 60000000 ~ "p", chr == "X" & end > 60000000 ~ "q",
    chr == "Y" & end < 11000000 ~ "p", chr == "Y" & end > 11000000 ~ "q",
  ), chr_arm = paste(chr, arm, sep = "."))

# Prepare data to calculate z-score per autosome
chrs <- paste(rep(as.character(1:22), each = 2), rep(c("p", "q"), 22), sep = ".")
sample_depths <- colSums(tpm)
chr_norm_expr <- matrix(0, length(chrs), ncol(tpm), 
                        dimnames = list(chrs, colnames(tpm)))
zscr <- chr_norm_expr

for(i in seq_along(chrs)){
  # Collect expression from each chromosome and normalise by sample depth
  chr_norm_expr[chrs[i], ] <- colSums(
    tpm[pull(filter(rna_genes, chr_arm == chrs[i]), Ensembl), ]
  ) / sample_depths
  # Compute z-score based on mean and sd across samples
  zscr[chrs[i], ] <- (chr_norm_expr[chrs[i], ] - mean(chr_norm_expr[chrs[i], ])) /
    sd(chr_norm_expr[chrs[i], ])
}

# Form the karyomap matrix using z-score > 1.65 for gain, < -1.65 for loss and 
# anything in-between for normal chromosome (z-score == 1.65 --> p-val ~ 0.05)
karyo <- zscr
karyo[zscr >= 1.65] <- 1
karyo[zscr <= -1.65] <- -1
karyo[zscr < 1.65 & zscr > -1.65] <- 0

# Create a heatmap
pheatmap(karyo, cluster_rows = FALSE, cluster_cols = FALSE)
