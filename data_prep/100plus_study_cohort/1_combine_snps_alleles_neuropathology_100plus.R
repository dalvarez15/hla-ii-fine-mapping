## =============================================================================
## Combine 100-plus Study SNP/HLA allele dosages with neuropathology and
## microglial marker data, linked by NBB and 100-plus Study ID
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   HUNDREDPLUS_NEUROPATHOLOGY_XLSX (raw_input_data/data_paths.R) - only the 100-plus/NBB
##     ID linkage is used from this file
##   HUNDREDPLUS_MICROGLIA_XLSX (raw_input_data/data_paths.R) - immunohistochemical
##     microglial marker quantification, Thal phase, Braak stage and CERAD score
##   data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv - combined HLA allele and SNP
##     dosages (from ../study_cohort_genetics/5_get_snps_alleles_dosages.R)
## Outputs:
##   data_prep/100plus_study_cohort/output/hla_alleles_snps_dosages_microglia_neuropathology_100plus.csv -
##     consumed by ../nbb_replication_cohort/3_microglia_neuropathology_regression.R
## =============================================================================

library(readxl)
library(data.table)

# Run from the repository root, e.g. `Rscript data_prep/100plus_study_cohort/1_combine_snps_alleles_neuropathology_100plus.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

dir.create("data_prep/100plus_study_cohort/output", recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

neuropathology <- read_excel(HUNDREDPLUS_NEUROPATHOLOGY_XLSX, skip = 1)
microglia_neuropathology <- read_excel(HUNDREDPLUS_MICROGLIA_XLSX)
alleles_snps_gwas <- fread("data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv", h = T, sep = ';', stringsAsFactors = F)

## -----------------------------------------------------------------------------
## 2. Link microglia, neuropathology and genetic data by sample ID
## -----------------------------------------------------------------------------

# Add the 100-plus Study ID ("ID") from the neuropathology overview into the
# microglia data (linked by NBB ID)
microglia_id <- merge(microglia_neuropathology, neuropathology[, c("NBB", "ID")], by = "NBB", all.x = TRUE)

# Add allele/SNP dosages (linked by NBB ID)
microglia_alleles_combined <- merge(microglia_id, alleles_snps_gwas, by.x = "NBB", by.y = "ID_NBB")

write.csv2(
    microglia_alleles_combined,
    "data_prep/100plus_study_cohort/output/hla_alleles_snps_dosages_microglia_neuropathology_100plus.csv",
    row.names = FALSE
)
