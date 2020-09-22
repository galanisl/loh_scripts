# Example execution: Rscript correct_amp_bias.R
# Note that this script corrects the VCF files produced by call_snps.sh in
# the vcf_filtered directory and outputs the corrected VCF files to the
# vcf_corrected folder, which should exist in the working directory.

library(vcfR)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# Get the fraction of reads supporting the Ref and Alt alleles in homozygous or
# heterozygous variants from a sample of interest
get_homhet_counts <- function(sample_file, var_type = "hom"){
  # Read VCF file
  vcf <- read.vcfR(sample_file)
  
  # Extract number of forward reads supporting the Ref and Alt alleles 
  fw_counts <- extract_gt_tidy(vcf, "ADF") %>% 
    separate(col = gt_ADF, into = c("ref_fw", "alt_fw"), sep = ",") %>% 
    mutate(ref_fw = as.integer(ref_fw), alt_fw = as.integer(alt_fw)) %>% 
    separate(col = gt_GT_alleles, into = c("ref", "alt"), sep = "/") %>% 
    select(Key, ref_fw, alt_fw, ref, alt)
  
  # Extract number of reverse reads supporting the Ref and Alt alleles   
  rv_counts <- extract_gt_tidy(vcf, "ADR") %>% 
    separate(col = gt_ADR, into = c("ref_rv", "alt_rv"), sep = ",") %>% 
    mutate(ref_rv = as.integer(ref_rv), alt_rv = as.integer(alt_rv)) %>% 
    select(Key, ref_rv, alt_rv)
  
  # Join the tables
  cts <- left_join(fw_counts, rv_counts, by = "Key") %>% 
    select(Key, ref_fw, ref_rv, alt_fw, alt_rv, ref, alt)
  
  # Focus on the 'var_type' calls and compute the fraction of reads supporting
  # the Ref and Alt alleles
  if(var_type == "hom"){
    cts <- cts %>% 
      filter(ref == alt)
  }else{
    cts <- cts %>% 
      filter(ref != alt)
  }
  cts <- cts %>% 
    mutate(tot = ref_fw + ref_rv + alt_fw + alt_rv,
           ref_tot = ref_fw + ref_rv, alt_tot = alt_fw + alt_rv,
           frac_ref = ref_tot/tot, frac_alt = alt_tot/tot,
           sample_name = sample_file) 
  
  return(cts)
}

# Transform hom-SNPs to het-SNPs based on a new threshold to call heterozygosity
make_het <- function(sample_file, threshold, output_folder = "vcf_corrected"){
  # Read VCF file
  vcf <- read.vcfR(sample_file)
  
  # Extract number of forward reads supporting the Ref and Alt alleles 
  fw_counts <- extract_gt_tidy(vcf, "ADF") %>% 
    separate(col = gt_ADF, into = c("ref_fw", "alt_fw"), sep = ",") %>% 
    mutate(ref_fw = as.integer(ref_fw), alt_fw = as.integer(alt_fw)) %>% 
    separate(col = gt_GT_alleles, into = c("ref", "alt"), sep = "/") %>% 
    select(Key, ref_fw, alt_fw, ref, alt)
  
  # Extract number of reverse reads supporting the Ref and Alt alleles   
  rv_counts <- extract_gt_tidy(vcf, "ADR") %>% 
    separate(col = gt_ADR, into = c("ref_rv", "alt_rv"), sep = ",") %>% 
    mutate(ref_rv = as.integer(ref_rv), alt_rv = as.integer(alt_rv)) %>% 
    select(Key, ref_rv, alt_rv)
  
  # Join the tables
  cts <- left_join(fw_counts, rv_counts, by = "Key") %>% 
    select(Key, ref_fw, ref_rv, alt_fw, alt_rv, ref, alt) %>% 
    mutate(tot = ref_fw + ref_rv + alt_fw + alt_rv,
           ref_tot = ref_fw + ref_rv, alt_tot = alt_fw + alt_rv,
           frac_ref = ref_tot/tot, frac_alt = alt_tot/tot,
           make_het_snp = ifelse(frac_ref >= threshold, TRUE, FALSE))
  
  # Correct
  if(sum(cts$make_het_snp) > 0){
  vcf@fix[cts$make_het_snp, "INFO"] <- str_replace(vcf@fix[cts$make_het_snp, "INFO"], 
                                                   "AC=2", "AC=1")
  }
  write.vcf(vcf, file = paste0(output_folder, "/", 
                               map_chr(strsplit(sample_file, "//"), `[`, 2), 
                               ".gz"))
}

vcf_files <- list.files("vcf_filtered/", full.names = TRUE)

# Correction of homozygous calls --------------------------------------------

hom_fractions <- map(vcf_files, get_homhet_counts, var_type = "hom")
names(hom_fractions) <- map_chr(strsplit(vcf_files, "//"), `[`, 2)

# Compute the median fraction of reads different from 0 that support the Ref 
# Allele in homozygous calls
med_fraction <- hom_fractions %>% 
  map_df(bind_rows) %>% 
  filter(frac_ref > 0) %>% 
  pull(frac_ref) %>% 
  median()

# Correct the homozygous calls based on the median fraction computed above
map_output <- map(vcf_files, make_het, threshold = med_fraction,
                  output_folder = "vcf_corrected/")

# Unzip the VCF files (region of interest is just 20kb so they're not that big)
system(paste0("gunzip ", "vcf_corrected/*.vcf.gz"))

