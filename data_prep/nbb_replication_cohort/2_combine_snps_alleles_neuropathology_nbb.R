## =============================================================================
## NBB replication cohort: HLA allele dosage matrix (two-field, DR-broad and
## one-field), merged with the haplotype lead SNP dosages and phenotypes from
## 1_snps_phenotypes.R. HLA-A/B/C and DPB1 were also imputed but are outside the
## scope of this HLA-II haplotype study and are not loaded here.
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   NBB_HLA_IMPUTATION_DIR (raw_input_data/data_paths.R) - HIBAG HLA imputation output per
##     locus (result_{DQA1,DQB1,DRB1}.txt), NBB cohort
##   data_prep/nbb_replication_cohort/output/nbb_snps_phenotypes.txt - QC'd haplotype lead
##     SNP dosages with phenotypes (from 1_snps_phenotypes.R)
## Outputs:
##   data_prep/nbb_replication_cohort/output/hla_imputation/imputation_probability_distribution.pdf -
##     QC plot: percentage of accepted genotype calls per locus
##   data_prep/nbb_replication_cohort/output/hla_imputation/hla_allele_dosages_nbb.csv -
##     NBB two-field + DR-broad + one-field allele dosage matrix
##   data_prep/nbb_replication_cohort/output/hla_alleles_snps_dosages_neuropathology_nbb.csv - combined SNP
##     + allele dosages + phenotypes, consumed by 3_microglia_neuropathology_regression.R
## =============================================================================

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

# Run from the repository root, e.g. `Rscript data_prep/nbb_replication_cohort/2_combine_snps_alleles_neuropathology_nbb.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

dir.create("data_prep/nbb_replication_cohort/output/hla_imputation", recursive = TRUE, showWarnings = FALSE)

# Posterior probability threshold for accepting an imputed genotype call,
# as recommended by HIBAG's authors (Zheng et al. 2014, Pharmacogenomics J.)
HIBAG_ACCEPT_PROB <- 0.5

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Reads and combines the per-locus HIBAG imputation result files in
# `data_path` into one long-format table, renaming the sample ID column.
read_hla_imputation <- function(data_path, id_col) {
    loci <- c("DQA1", "DQB1", "DRB1")
    imputed_hla <- do.call(rbind, lapply(loci, function(locus) {
        d <- read.delim(paste0(data_path, "/result_", locus, ".txt"))
        d$locus <- locus
        d
    }))
    colnames(imputed_hla)[colnames(imputed_hla) == "sample.id"] <- id_col
    imputed_hla
}

# Adds a grouped allele-count column (sum of dosages across every column
# matching `pattern`, e.g. every DRB1 two-field allele within one ancestral
# structure), NA where none of the matching alleles were called.
add_grouped_allele_column <- function(dt, new_col, pattern) {
    dt[, (new_col) := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns(pattern)]
    dt
}

## -----------------------------------------------------------------------------
## 1. Load and combine per-locus imputation results
## -----------------------------------------------------------------------------

imputation_nbb <- read_hla_imputation(NBB_HLA_IMPUTATION_DIR, "ID_NBB")
imputation_nbb$ID_NBB <- paste0("NBB ", sapply(strsplit(imputation_nbb$ID_NBB, "_"), "[", 1))

## -----------------------------------------------------------------------------
## 2. QC plot: accepted genotype calls per locus
## -----------------------------------------------------------------------------

qc_plot <- imputation_nbb %>%
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
ggsave("data_prep/nbb_replication_cohort/output/hla_imputation/imputation_probability_distribution.pdf", plot = qc_plot, width = 10, height = 6)

## -----------------------------------------------------------------------------
## 3. Build the two-field allele dosage matrix
## -----------------------------------------------------------------------------

allele_counts <- imputation_nbb %>%
    filter(prob > HIBAG_ACCEPT_PROB) %>%
    select(ID_NBB, locus, allele1, allele2) %>%
    pivot_longer(cols = c(allele1, allele2), names_to = "allele_num", values_to = "allele_code") %>%
    mutate(allele = paste0(locus, "*", allele_code)) %>%
    group_by(ID_NBB, locus, allele) %>%
    summarise(count = n(), .groups = "drop")

all_alleles <- imputation_nbb %>%
    select(locus, allele1, allele2) %>%
    pivot_longer(cols = c(allele1, allele2), values_to = "allele_code") %>%
    mutate(allele = paste0(locus, "*", allele_code)) %>%
    distinct(locus, allele)

all_combinations <- expand.grid(
    ID_NBB = unique(imputation_nbb$ID_NBB),
    locus = unique(imputation_nbb$locus),
    allele = unique(all_alleles$allele),
    stringsAsFactors = FALSE
) %>%
    inner_join(all_alleles, by = c("locus", "allele"))

allele_counts_full <- all_combinations %>%
    left_join(allele_counts, by = c("ID_NBB", "locus", "allele")) %>%
    left_join(
        imputation_nbb %>%
            group_by(ID_NBB, locus) %>%
            summarise(accepted = all(prob > HIBAG_ACCEPT_PROB), .groups = "drop"),
        by = c("ID_NBB", "locus")
    ) %>%
    mutate(count = ifelse(!accepted, NA, replace_na(count, 0))) %>%
    select(ID_NBB, locus, allele, count)

allele_count_matrix <- allele_counts_full %>%
    select(ID_NBB, allele, count) %>%
    pivot_wider(names_from = allele, values_from = count)
setDT(allele_count_matrix)

## -----------------------------------------------------------------------------
## 4. Derive DR-broad ancestral structures and one-field allele columns
## -----------------------------------------------------------------------------

add_grouped_allele_column(allele_count_matrix, "DR1", "DRB1\\*01|DRB1\\*10")
add_grouped_allele_column(allele_count_matrix, "DR8", "DRB1\\*08")           # derivative of DR3
add_grouped_allele_column(allele_count_matrix, "DR2", "DRB1\\*15|DRB1\\*16")
add_grouped_allele_column(allele_count_matrix, "DR3", "DRB1\\*03|DRB1\\*11|DRB1\\*12|DRB1\\*13|DRB1\\*14")
add_grouped_allele_column(allele_count_matrix, "DR4", "DRB1\\*04|DRB1\\*07|DRB1\\*09")
add_grouped_allele_column(allele_count_matrix, "DRB1*04", "DRB1\\*04") # one-field DRB1*04, kept separately from the DR4 structure above

haplotype_cols <- c("DR1", "DR8", "DR2", "DR3", "DR4", "DRB1*04")
allele_count_matrix[, (haplotype_cols) := lapply(.SD, as.integer), .SDcols = haplotype_cols]
stopifnot(all(sapply(haplotype_cols, function(col) all(allele_count_matrix[[col]] %in% c(0L, 1L, 2L, NA)))))

for (f in c("01", "10", "08", "15", "16", "03", "11", "12", "13", "14", "07", "09")) {
    add_grouped_allele_column(allele_count_matrix, paste0("DRB1*", f), paste0("DRB1\\*", f))
}

for (g in c("DQB1", "DQA1")) {
    g_cols <- grep(paste0("^", g, "\\*"), colnames(allele_count_matrix), value = TRUE)
    g_first <- unique(sub(paste0("^", g, "\\*([0-9]+).*"), "\\1", g_cols))
    for (f in g_first) {
        new_col <- paste0(g, "*", f)
        patt <- paste0(g, "\\*", f, "(:|$)") # ensure distinct matching
        add_grouped_allele_column(allele_count_matrix, new_col, patt)
        allele_count_matrix[, (new_col) := as.integer(get(new_col))]
    }
}

write.csv(
    allele_count_matrix,
    "data_prep/nbb_replication_cohort/output/hla_imputation/hla_allele_dosages_nbb.csv",
    row.names = FALSE
)

## -----------------------------------------------------------------------------
## 5. Merge with haplotype lead SNP dosages and phenotypes
## -----------------------------------------------------------------------------

nbb_snps_phenotypes <- read.table(
    "data_prep/nbb_replication_cohort/output/nbb_snps_phenotypes.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

merged_nbb_alleles_snps <- merge(
    nbb_snps_phenotypes,
    as.data.frame(allele_count_matrix),
    by.x = "IID", by.y = "ID_NBB"
)

# HLA imputation QC, restricted to this final cohort (cited in the manuscript)
imputation_final_cohort <- imputation_nbb[imputation_nbb$ID_NBB %in% merged_nbb_alleles_snps$IID, ]
n_accepted <- sum(imputation_final_cohort$prob > HIBAG_ACCEPT_PROB)
n_total <- nrow(imputation_final_cohort)
message(sprintf(
    "HLA imputation QC (NBB): %d of %d (%.1f%%) imputed HLA genotypes accepted (posterior probability > %.1f) across HLA-DRB1, HLA-DQA1 and HLA-DQB1 in %d individuals",
    n_accepted, n_total, 100 * n_accepted / n_total, HIBAG_ACCEPT_PROB, length(unique(imputation_final_cohort$ID_NBB))
))

write.csv2(
    merged_nbb_alleles_snps,
    "data_prep/nbb_replication_cohort/output/hla_alleles_snps_dosages_neuropathology_nbb.csv",
    row.names = FALSE
)
