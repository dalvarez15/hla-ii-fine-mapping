## =============================================================================
## Combine per-locus HLA imputation results (DQA1, DQB1, DRB1) into an allele
## dosage matrix. HLA-A/B/C and DPB1 were also imputed but are outside the
## scope of this HLA-II haplotype study and are not loaded here.
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   STUDY_COHORT_HLA_IMPUTATION_DIR (raw_input_data/data_paths.R) - HIBAG HLA imputation
##     output per locus (result_{DQA1,DQB1,DRB1}.txt), study cohort
## Outputs:
##   data_prep/study_cohort_genetics/output/hla_imputation/imputation_probability_distribution.pdf -
##     QC plot: percentage of accepted genotype calls per locus
##   data_prep/study_cohort_genetics/output/hla_imputation/hla_allele_dosages_study_cohort.csv -
##     allele dosage matrix, consumed by 5_get_snps_alleles_dosages.R
## =============================================================================

library(ggplot2)
library(tidyr)
library(dplyr)

# Run from the repository root, e.g. `Rscript data_prep/study_cohort_genetics/4_read_hla_imputation.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

OUTDIR <- "data_prep/study_cohort_genetics/output/hla_imputation"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Posterior probability threshold for accepting an imputed genotype call,
# as recommended by HIBAG's authors (Zheng et al. 2014, Pharmacogenomics J.)
HIBAG_ACCEPT_PROB <- 0.5

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Reads and combines the per-locus HIBAG imputation result files in
# `data_path` into one long-format table (one row per participant x locus).
read_hla_imputation <- function(data_path) {
    loci <- c("DQA1", "DQB1", "DRB1")
    imputed_hla <- do.call(rbind, lapply(loci, function(locus) {
        d <- read.delim(paste0(data_path, "/result_", locus, ".txt"))
        d$locus <- locus
        d
    }))
    colnames(imputed_hla)[colnames(imputed_hla) == "sample.id"] <- "ID_GWAS"
    imputed_hla
}

## -----------------------------------------------------------------------------
## 1. Load and combine per-locus imputation results
## -----------------------------------------------------------------------------

imputation_combined <- read_hla_imputation(STUDY_COHORT_HLA_IMPUTATION_DIR)

## -----------------------------------------------------------------------------
## 2. QC plot: accepted genotype calls per locus
## -----------------------------------------------------------------------------

qc_plot <- imputation_combined %>%
    group_by(locus) %>%
    summarise(
        percent_accepted = round(100 * sum(prob > HIBAG_ACCEPT_PROB) / n(), 1)
    ) %>%
    ggplot(aes(x = locus, y = percent_accepted, fill = locus)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = paste0(percent_accepted, "%")), vjust = -0.5, size = 5) +
    labs(title = "Percentage of Accepted Imputation Rows per Locus (prob > HIBAG_ACCEPT_PROB)",
         x = "Locus",
         y = "Accepted (%)") +
    theme_minimal() +
    ylim(0, 100)
ggsave(file.path(OUTDIR, "imputation_probability_distribution.pdf"), plot = qc_plot, width = 10, height = 6)

## -----------------------------------------------------------------------------
## 3. Build the allele dosage matrix
## -----------------------------------------------------------------------------

# Count of each two-field allele per participant, among accepted calls
allele_counts <- imputation_combined %>%
    filter(prob > HIBAG_ACCEPT_PROB) %>%
    select(ID_GWAS, locus, allele1, allele2) %>%
    pivot_longer(cols = c(allele1, allele2), names_to = "allele_num", values_to = "allele_code") %>%
    mutate(allele = paste0(locus, "*", allele_code)) %>%
    group_by(ID_GWAS, locus, allele) %>%
    summarise(count = n(), .groups = "drop")

# Every possible allele per locus
all_alleles <- imputation_combined %>%
    select(locus, allele1, allele2) %>%
    pivot_longer(cols = c(allele1, allele2), values_to = "allele_code") %>%
    mutate(allele = paste0(locus, "*", allele_code)) %>%
    distinct(locus, allele)

# Full grid of participant x locus x allele, so alleles never carried by a
# participant appear as an explicit 0 (rather than being silently dropped),
# while loci that failed the probability threshold stay NA
all_combinations <- expand.grid(
    ID_GWAS = unique(imputation_combined$ID_GWAS),
    locus = unique(imputation_combined$locus),
    allele = unique(all_alleles$allele),
    stringsAsFactors = FALSE
) %>%
    inner_join(all_alleles, by = c("locus", "allele"))

allele_counts_full <- all_combinations %>%
    left_join(allele_counts, by = c("ID_GWAS", "locus", "allele")) %>%
    left_join(
        imputation_combined %>%
            group_by(ID_GWAS, locus) %>%
            summarise(accepted = all(prob > HIBAG_ACCEPT_PROB), .groups = "drop"),
        by = c("ID_GWAS", "locus")
    ) %>%
    mutate(count = ifelse(!accepted, NA, replace_na(count, 0))) %>%
    select(ID_GWAS, locus, allele, count)

allele_count_matrix <- allele_counts_full %>%
    select(ID_GWAS, allele, count) %>%
    pivot_wider(names_from = allele, values_from = count)

write.csv(
    allele_count_matrix,
    file.path(OUTDIR, "hla_allele_dosages_study_cohort.csv"),
    row.names = FALSE
)
