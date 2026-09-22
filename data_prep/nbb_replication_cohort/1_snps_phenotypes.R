## =============================================================================
## NBB replication cohort: dosages of the three haplotype lead SNPs carried
## forward for replication (Hap-R, Hap-B, Hap-Y), re-oriented to their minor
## allele and verified against the primary cohort's COJO haplotype table,
## merged with NBB phenotypes and restricted to participants passing genetic
## QC (European ancestry, unrelated, age/sex/post-mortem-delay available)
## with an AD or control neuropathological diagnosis.
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   NBB_CHR6_DOSAGE_RAW (raw_input_data/data_paths.R) - PLINK dosage export for chr6, NBB
##     genotyping array
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - COJO haplotype
##     lead SNPs, minor-allele oriented (from ../study_cohort_genetics/3_cojo_annotate.R)
##   NBB_PHENOTYPES_RDATA (raw_input_data/data_paths.R) - NBB sample phenotypes and QC flags
## Outputs:
##   data_prep/nbb_replication_cohort/output/nbb_snps_phenotypes.txt - QC'd haplotype lead SNP
##     dosages with phenotypes, consumed by 2_combine_snps_alleles_neuropathology_nbb.R
## =============================================================================

library(dplyr)

# Run from the repository root, e.g. `Rscript data_prep/nbb_replication_cohort/1_snps_phenotypes.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

dir.create("data_prep/nbb_replication_cohort/output", recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. Load SNP dosages and restrict to the three replication haplotype lead SNPs
## -----------------------------------------------------------------------------

snps <- read.table(NBB_CHR6_DOSAGE_RAW, h = T, check.names = FALSE, stringsAsFactors = F)
colnames(snps)[7:ncol(snps)] <- paste0("chr", colnames(snps)[7:ncol(snps)]) # add chr prefix to match cojo_ann

cojo_ann <- read.csv2("data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv", header = TRUE, check.names = FALSE)
cojo_ann_paper <- cojo_ann[cojo_ann$snp_easy %in% c("Hap-R", "Hap-B", "Hap-Y"), ]

snps_cojo <- snps[, c("IID", unlist(lapply(cojo_ann_paper$SNP, function(snp) grep(snp, colnames(snps), value = TRUE)))), drop = FALSE]

## -----------------------------------------------------------------------------
## 2. Re-orient to the minor allele and verify against the COJO table
## -----------------------------------------------------------------------------

# In the NBB array export, all three SNPs are coded on the major allele;
# flip both the column names' trailing allele suffix and the dosages
colnames(snps_cojo)[which(colnames(snps_cojo) == "chr6:32447376:C:T_C")] <- "chr6:32447376:C:T_T"
colnames(snps_cojo)[which(colnames(snps_cojo) == "chr6:32592593:G:T_G")] <- "chr6:32592593:G:T_T"
colnames(snps_cojo)[which(colnames(snps_cojo) == "chr6:32609869:A:G_G")] <- "chr6:32609869:A:G_A"
snps_cojo[, -1] <- lapply(snps_cojo[, -1, drop = FALSE], function(x) 2 - x)

calc_freq <- colSums(snps_cojo[, -1, drop = FALSE], na.rm = TRUE) / (2 * nrow(snps_cojo))
freq_mismatch <- abs(calc_freq[match(cojo_ann_paper$snp_dos, names(calc_freq))] - cojo_ann_paper$freq_geno) > 0.05
if (any(freq_mismatch, na.rm = TRUE)) {
    stop("NBB allele frequency does not match the primary cohort's COJO table for: ", paste(cojo_ann_paper$snp_dos[freq_mismatch], collapse = ", "))
}

## -----------------------------------------------------------------------------
## 3. Merge phenotypes and restrict to QC-passing AD/control samples
## -----------------------------------------------------------------------------

load(NBB_PHENOTYPES_RDATA)
nbb_phenotypes <- dt
snps_cojo$IID <- paste0("NBB ", sapply(strsplit(snps_cojo$IID, "_"), "[", 1))
nbb_snps_phenotypes <- merge(snps_cojo, nbb_phenotypes, by = "IID")

nbb_snps_phenotypes <- nbb_snps_phenotypes %>%
  filter(
    european == TRUE,           # European ancestry only
    is_unrelated == TRUE,       # exclude related individuals
    age != -1,                  # age available
    most_likely_sex != "Unknown",
    !is.na(post_mortem_delay),
    nd_current_mapped %in% c("AD", "CON") # AD or control neuropathological diagnosis
  )

write.table(nbb_snps_phenotypes, "data_prep/nbb_replication_cohort/output/nbb_snps_phenotypes.txt", row.names = FALSE, quote = FALSE, sep = "\t")
