## =============================================================================
## Hap-B/Hap-R/Hap-Y vs. microglial markers and neuropathology, 100-plus Study
## and NBB replication cohort
##
## For every genetic feature of Hap-B, Hap-R and Hap-Y, and every
## microglial-marker/neuropathology outcome, fits a linear regression
## (100-plus: correcting for genetic PCs 1-5, sex and age at death; NBB:
## additionally correcting for post-mortem delay) and applies an
## effective-number-of-tests (Li & Ji) FDR correction within each outcome
## class. This is the one data_prep/ step spanning two cohorts, so it runs
## last, after both ../100plus_study_cohort/1_combine_snps_alleles_neuropathology_100plus.R
## and ../nbb_replication_cohort/2_combine_snps_alleles_neuropathology_nbb.R.
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - COJO haplotype
##     lead SNPs (minor-allele oriented), mapping dosage columns to easy-to-read labels; from
##     ../study_cohort_genetics/3_cojo_annotate.R
##   data_prep/100plus_study_cohort/output/hla_alleles_snps_dosages_microglia_neuropathology_100plus.csv - combined
##     microglial-marker, neuropathology, HLA allele/SNP and covariate data for the 100-plus
##     Study post-mortem cohort; from ../100plus_study_cohort/1_combine_snps_alleles_neuropathology_100plus.R
##   data_prep/nbb_replication_cohort/output/hla_alleles_snps_dosages_neuropathology_nbb.csv - NBB HLA
##     allele/SNP and neuropathology phenotype data; from 2_combine_snps_alleles_neuropathology_nbb.R
##   NBB_PCS_EIGENVEC (raw_input_data/data_paths.R) - genetic PCs for NBB samples; no
##     generating script in this repository
## Outputs:
##   data_prep/nbb_replication_cohort/output/microglia_neuropathology_regression.csv - full
##     association statistics (100-plus + NBB), consumed by table3_supptable2.R
## =============================================================================

library(data.table)
library(dplyr)

# Run from the repository root, e.g. `Rscript data_prep/nbb_replication_cohort/3_microglia_neuropathology_regression.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Renames columns of `df` according to a named vector `mapping`
# (names = current column names, values = new names); columns not present
# in `df` are silently skipped.
rename_columns <- function(df, mapping) {
  hits <- intersect(names(mapping), colnames(df))
  colnames(df)[match(hits, colnames(df))] <- mapping[hits]
  df
}

# Adds an easy-to-read SNP column (snp_easy, e.g. "rs35472547") to `df` for
# every dosage SNP column (snp_dos) present in `df`, based on the COJO
# haplotype lookup table `snp_lookup`.
add_easy_snp_labels <- function(df, snp_lookup) {
  for (i in seq_len(nrow(snp_lookup))) {
    dos_col  <- snp_lookup$snp_dos[i]
    easy_col <- snp_lookup$snp_easy[i]
    if (!is.null(dos_col) && dos_col %in% colnames(df) && !(easy_col %in% colnames(df))) {
      df[[easy_col]] <- df[[dos_col]]
    }
  }
  df
}

# Fits one linear regression per (outcome, genetic feature) pair:
#   outcome ~ feature + covariates
# and extracts the feature's beta, SE, standardized beta, p-value and N.
#
# `outcome_labels` optionally relabels each outcome in the output (e.g. to
# note which other outcome it was additionally corrected for).
run_hla_regression_loop <- function(data, outcomes, features, formula_covariates,
                                     outcome_labels = setNames(outcomes, outcomes)) {
  results <- data.frame()
  for (test in outcomes) {
    for (feature in features) {
      dummy <- paste0("snp_dummy_", make.names(feature))
      data[[dummy]] <- data[[feature]]

      model_vars <- c(test, dummy, formula_covariates)
      idx <- complete.cases(data[, model_vars])

      if (sum(idx) > 10) {
        df_model <- data[idx, ]

        rhs <- paste(paste0("`", c(dummy, formula_covariates), "`"), collapse = " + ")
        fml <- as.formula(paste0("`", test, "` ~ ", rhs))

        fit <- lm(fml, data = df_model)
        coef_summary <- summary(fit)$coefficients

        if (dummy %in% rownames(coef_summary)) {
          p_col <- grep("^Pr\\(", colnames(coef_summary), value = TRUE)[1]

          beta <- coef_summary[dummy, "Estimate"]
          pval <- coef_summary[dummy, p_col]

          sd_x <- sd(df_model[[dummy]], na.rm = TRUE)
          sd_y <- sd(df_model[[test]], na.rm = TRUE)
          beta_std <- ifelse(sd_x == 0 | sd_y == 0, NA, beta * (sd_x / sd_y))

          results <- rbind(results, data.frame(
            test_var  = outcome_labels[[test]],
            snp_var   = feature,
            beta      = beta,
            std_error = coef_summary[dummy, "Std. Error"],
            beta_std  = beta_std,
            pval      = pval,
            n         = sum(idx)
          ))
        }
      }
      data[[dummy]] <- NULL
    }
  }
  results
}

# Effective number of independent tests among correlated markers (Li & Ji,
# 2005): sum of the eigenvalues of the marker correlation matrix, capped at 1
# each, used below to correct FDR for the correlation between the haplotypes'
# SNP/allele features.
calc_meff_li_ji <- function(X) {
  R <- cor(X, use = "pairwise.complete.obs")
  ev <- eigen(R, only.values = TRUE)$values
  sum(pmin(ev, 1))
}

# FDR-corrects p-values within each outcome class (e.g. microglial markers
# vs. neuropathology scores), scaling the Benjamini-Hochberg correction by
# the effective (Li & Ji) rather than nominal number of markers tested.
add_meff_fdr <- function(results, outcome_groups, m_eff_per_marker) {
  results$p_fdr <- NA_real_
  for (group in outcome_groups) {
    idx <- results$test_var %in% group
    pvals <- results$pval[idx]

    n_pheno <- length(unique(results$test_var[idx]))
    m_eff_total <- m_eff_per_marker * n_pheno

    q_bh <- p.adjust(pvals, method = "fdr")
    q_meff <- pmin(q_bh * (m_eff_total / length(pvals)), 1)

    results$p_fdr[idx] <- q_meff
  }
  results
}

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

microglia_alleles_combined <- read.csv2(
  "data_prep/100plus_study_cohort/output/hla_alleles_snps_dosages_microglia_neuropathology_100plus.csv",
  check.names = FALSE
)

cojo_haplotypes <- read.csv2(
  "data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv"
)

microglia_alleles_combined <- add_easy_snp_labels(microglia_alleles_combined, cojo_haplotypes)

# Display labels: lead SNPs by rs number, and outcome columns renamed for
# readability (Iba1, HLA class II, Greek letters for amyloid-beta)
outcome_and_feature_display_labels <- c(
  "Hap-R"               = "rs9469112",
  "Hap-B"               = "rs35472547",
  "Hap-Y"               = "rs4335021",
  "Iba-1  load"         = "Iba1 load",
  "HLA-DR load"         = "HLA class II load",
  "Braak stage"         = "Braak NFT stage (100+)",
  "Thal AB stage"       = "Thal Aβ phase",
  "AB load"             = "Aβ load",
  "AB40 load"           = "Aβ40 load",
  "AB42 load"           = "Aβ42 load"
)
microglia_alleles_combined <- rename_columns(microglia_alleles_combined, outcome_and_feature_display_labels)

## -----------------------------------------------------------------------------
## 2. Define outcomes and haplotype genetic features (100-plus cohort)
## -----------------------------------------------------------------------------

neuropathology_scores <- c("CERAD score", "Braak NFT stage (100+)", "Thal Aβ phase")
neuropathology_loads  <- c("Aβ load", "Aβ40 load", "Aβ42 load", "pTau217 load", "AT8 load", "GT-38 load")
microglia_loads       <- c("HLA class II load", "CD11c load", "CD68 load", "Iba1 load", "P2RY12 load")
neuropathology_all    <- c(neuropathology_scores, neuropathology_loads)

# Genetic features of each haplotype: lead SNP (rs), DR-broad code, and
# two-field HLA-DRB1/DQA1/DQB1 alleles
haplotype_features <- c(
  "rs35472547", "DR4", "DRB1*04:01", "DQA1*03:01", "DQB1*03:02", # Hap-B (blue)
  "rs9469112",  "DR1", "DRB1*01:01", "DQA1*01:01", "DQB1*05:01", # Hap-R (red)
  "rs4335021",  "DR2", "DRB1*15:01", "DQA1*01:02", "DQB1*06:02"  # Hap-Y (yellow)
)

haplotype_lookup <- c(
  rs35472547 = "Hap-B", DR4 = "Hap-B", `DRB1*04:01` = "Hap-B", `DQA1*03:01` = "Hap-B", `DQB1*03:02` = "Hap-B",
  rs9469112 = "Hap-R", DR1 = "Hap-R", `DRB1*01:01` = "Hap-R", `DQA1*01:01` = "Hap-R", `DQB1*05:01` = "Hap-R",
  rs4335021 = "Hap-Y", DR2 = "Hap-Y", `DRB1*15:01` = "Hap-Y", `DQA1*01:02` = "Hap-Y", `DQB1*06:02` = "Hap-Y"
)

# 100-plus cohort covariates. Post-mortem delay is not available for the AD
# cases, so it is left out here rather than corrected for, to keep all 96
# participants (89 centenarians + 7 AD cases) in every model.
pc_covariates <- paste0("PC", 1:5)
cohort100plus_formula_covariates <- c(pc_covariates, "Sex", "Age_at_death")

## -----------------------------------------------------------------------------
## 3. Regression models: 100-plus microglia and neuropathology
## -----------------------------------------------------------------------------

results <- run_hla_regression_loop(
  data = microglia_alleles_combined,
  outcomes = c(neuropathology_all, microglia_loads),
  features = haplotype_features,
  formula_covariates = cohort100plus_formula_covariates
)
results$snp_var <- factor(results$snp_var, levels = haplotype_features)

# Effective number of independent tests among the haplotype features, used
# for the Li & Ji-scaled FDR correction below
meff_100plus <- calc_meff_li_ji(microglia_alleles_combined[, haplotype_features])

results <- add_meff_fdr(results, list(neuropathology_all, microglia_loads), meff_100plus)
results$test_var <- factor(results$test_var, levels = c(microglia_loads, neuropathology_all))

## -----------------------------------------------------------------------------
## 4. NBB replication cohort: sample overlap and regression models
## -----------------------------------------------------------------------------

merged_nbb_alleles_snps <- read.csv2(
  "data_prep/nbb_replication_cohort/output/hla_alleles_snps_dosages_neuropathology_nbb.csv",
  check.names = FALSE
)
merged_nbb_alleles_snps <- add_easy_snp_labels(merged_nbb_alleles_snps, cojo_haplotypes)

# Display labels for NBB lead SNPs; DR-broad major MHC-II structures are
# labeled DR1/DR2/DR4 directly in the input.
merged_nbb_alleles_snps <- rename_columns(merged_nbb_alleles_snps, c(
  "Hap-R" = "rs9469112",
  "Hap-B" = "rs35472547",
  "Hap-Y" = "rs4335021"
))

# Convert Braak amyloid-beta stage from letter (O/A/B/C) to numeric (0-3)
braak_map <- c("O" = 0, "A" = 1, "B" = 2, "C" = 3)
merged_nbb_alleles_snps$braak_ab <- braak_map[as.character(merged_nbb_alleles_snps$braak_ab)]
merged_nbb_alleles_snps <- rename_columns(merged_nbb_alleles_snps, c(
  "braak_ab"   = "Braak amyloid-β stage",
  "braak_nft"  = "Braak NFT stage (NBB)",
  "braak_lewy" = "Braak α-synuclein stage"
))
nbb_outcomes <- c("Braak amyloid-β stage", "Braak NFT stage (NBB)", "Braak α-synuclein stage")

pcs <- read.table(NBB_PCS_EIGENVEC, col.names = c("FID", "IID", paste0("PC", 1:10)), stringsAsFactors = FALSE)
pcs$IID <- paste0("NBB ", pcs$IID)
merged_nbb_alleles_snps_test <- merge(merged_nbb_alleles_snps, pcs, by = "IID")

nbb_covariates <- c("most_likely_sex", "age", "post_mortem_delay", pc_covariates)

# Main NBB regressions: each of the three neuropathology outcomes against
# every haplotype feature
nbb_results <- run_hla_regression_loop(
  data = merged_nbb_alleles_snps_test,
  outcomes = nbb_outcomes,
  features = haplotype_features,
  formula_covariates = nbb_covariates
)

# Braak amyloid-beta and Braak NFT stage are additionally cross-adjusted for
# one another, since amyloid and tau pathology are correlated
nbb_outcomes_tau_amyloid <- c("Braak amyloid-β stage", "Braak NFT stage (NBB)")
nbb_results_tau_amyloid <- data.frame()
for (test in nbb_outcomes_tau_amyloid) {
  other_outcome <- setdiff(nbb_outcomes_tau_amyloid, test)
  nbb_results_tau_amyloid <- rbind(
    nbb_results_tau_amyloid,
    run_hla_regression_loop(
      data = merged_nbb_alleles_snps_test,
      outcomes = test,
      features = haplotype_features,
      formula_covariates = c(nbb_covariates, other_outcome),
      outcome_labels = setNames(paste0(test, "\ncorrected for:\n", other_outcome), test)
    )
  )
}
nbb_results <- rbind(nbb_results, nbb_results_tau_amyloid)
nbb_results$snp_var <- factor(nbb_results$snp_var, levels = haplotype_features)

meff_nbb <- calc_meff_li_ji(merged_nbb_alleles_snps_test[, haplotype_features])
nbb_results <- add_meff_fdr(nbb_results, list(unique(nbb_results$test_var)), meff_nbb)
nbb_results$test_var <- factor(nbb_results$test_var, levels = unique(nbb_results$test_var))

## -----------------------------------------------------------------------------
## 5. Combine and write the full association statistics
## -----------------------------------------------------------------------------

format_association_table <- function(regression_results) {
  regression_results %>%
    mutate(
      Haplotype = unname(haplotype_lookup[as.character(snp_var)]),
      `Outcome variable` = test_var,
      Beta = beta,
      SE = std_error,
      `Std. Beta` = beta_std,
      P = pval,
      `FDR P` = p_fdr,
      N = n
    ) %>%
    select(Haplotype, Feature = snp_var, `Outcome variable`, Beta, SE, `Std. Beta`, P, `FDR P`, N)
}

association_stats <- bind_rows(
  format_association_table(results),
  format_association_table(nbb_results)
) %>%
  arrange(factor(Haplotype, levels = c("Hap-B", "Hap-R", "Hap-Y")), Feature, `Outcome variable`)
association_stats$`Outcome variable` <- gsub("\n", " ", association_stats$`Outcome variable`)

fwrite(
  association_stats,
  file = "data_prep/nbb_replication_cohort/output/microglia_neuropathology_regression.csv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
