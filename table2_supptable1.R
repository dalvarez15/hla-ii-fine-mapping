## =============================================================================
## Supplementary Table 1 & Table 2: association of HLA-II haplotype genetic
## features (lead SNP, major MHC-II structure, one-field and two-field HLA
## alleles) with CHC/control/AD status, for Hap-B, Hap-R and Hap-Y. Table 2 is
## the subset reaching nominal significance (P<0.05) in at least one of the
## three pairwise comparisons (CHC vs. control, CHC vs. AD, control vs. AD).
## =============================================================================
##
## Inputs:
##   data_prep/study_cohort_genetics/output/hla_alleles_snps_regression.csv - regression
##     results for HLA-II haplotype genetic features: allele frequencies, sample sizes, and
##     CHC/control/AD association statistics; from 7_allele_snp_regression.R
## Outputs:
##   tables/supptable1.tsv
##   tables/supptable1_rounded.tsv
##   tables/table2.tsv
## =============================================================================

library(data.table)
library(dplyr)
library(tidyr)

# Run from the repository root, e.g. `Rscript table2_supptable1.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Builds one row per genetic feature of a haplotype, with sample size, carrier
# count/frequency in each group (CHC, control, AD), and association statistics
# (beta, SE, P, OR with 95% CI) for each of the three pairwise comparisons.
#
# Association labels differ slightly between resolution levels in the input
# file (e.g. "cent_cntr" vs. "cent_ctrl"), so the comparison of interest is
# resolved from a list of candidate labels rather than assumed fixed.
#
# Input:
#   regression_df - long-format regression results (one row per
#                    allele x association), with columns allele, association,
#                    beta, standard_error, p_value, odds_ratio, lower_ci,
#                    upper_ci, num_cases, num_controls, count_cases,
#                    count_controls, freq_cases, freq_controls
#   hap_name       - display label for the haplotype (e.g. "Hap-B")
#   features       - character vector of allele/SNP labels belonging to this
#                    haplotype, at every resolution level
# Output: one row per feature, with columns Haplotype, Feature, N.*, Count.*,
#   Freq.* (per group) and Beta/SE/P/ORCI per pairwise comparison.
build_haplotype_table <- function(regression_df, hap_name, features) {
  available_assocs <- unique(regression_df$association)

  resolve_assoc <- function(candidates, label) {
    hit <- intersect(candidates, available_assocs)
    if (length(hit) == 0) {
      stop(paste0("Missing association for ", label, ". Tried: ", paste(candidates, collapse = ", ")))
    }
    hit[1]
  }

  assoc_chc_ad <- resolve_assoc(c("cent_ad"), "CHC_AD")
  assoc_chc_cntr <- resolve_assoc(c("cent_cntr", "cent_ctrl", "centr_cntr"), "CHC_CNTR")
  assoc_cntr_ad <- resolve_assoc(c("ad_cntr", "ctrl_ad"), "CNTR_AD")

  regs <- regression_df %>%
    filter(allele %in% features,
           association %in% c(assoc_chc_ad, assoc_chc_cntr, assoc_cntr_ad)) %>%
    select(
      allele, association,
      beta, standard_error, p_value, odds_ratio, lower_ci, upper_ci,
      num_cases, num_controls,
      count_cases, count_controls,
      freq_cases, freq_controls
    )

  base_df <- data.frame(Feature = features, stringsAsFactors = FALSE)
  if (nrow(regs) == 0) {
    return(base_df %>% mutate(Haplotype = hap_name))
  }

  wide <- regs %>%
    pivot_wider(
      names_from = association,
      values_from = c(
        beta, standard_error, p_value, odds_ratio, lower_ci, upper_ci,
        num_cases, num_controls,
        count_cases, count_controls,
        freq_cases, freq_controls
      ),
      names_sep = "_"
    )

  joined <- base_df %>%
    left_join(wide, by = c("Feature" = "allele"))

  safe_get <- function(df, col_name) {
    if (col_name %in% names(df)) {
      return(df[[col_name]])
    }
    rep(NA_real_, nrow(df))
  }

  fmt_orci <- function(or, lci, uci) {
    ifelse(
      is.na(or) | is.na(lci) | is.na(uci),
      NA_character_,
      paste0(sprintf("%.2f", or), " [", sprintf("%.2f", lci), "-", sprintf("%.2f", uci), "]")
    )
  }

  d <- joined

  # Group labels follow the two comparisons that touch each group: e.g. the
  # CHC count/frequency is reported by the CHC-vs-control model where CHC is
  # coded as "controls" in that model, falling back to the CHC-vs-AD model
  # (where CHC is coded as "controls") when the former is missing.
  d %>%
    transmute(
      Haplotype = hap_name,
      Feature = Feature,

      N.CHC = coalesce(
        safe_get(d, paste0("num_controls_", assoc_chc_ad)),
        safe_get(d, paste0("num_controls_", assoc_chc_cntr))
      ),
      N.CNTR = coalesce(
        safe_get(d, paste0("num_cases_", assoc_chc_cntr)),
        safe_get(d, paste0("num_controls_", assoc_cntr_ad))
      ),
      N.AD = coalesce(
        safe_get(d, paste0("num_cases_", assoc_chc_ad)),
        safe_get(d, paste0("num_cases_", assoc_cntr_ad))
      ),

      Count.CHC = coalesce(
        safe_get(d, paste0("count_controls_", assoc_chc_ad)),
        safe_get(d, paste0("count_controls_", assoc_chc_cntr))
      ),
      Count.CNTR = coalesce(
        safe_get(d, paste0("count_cases_", assoc_chc_cntr)),
        safe_get(d, paste0("count_controls_", assoc_cntr_ad))
      ),
      Count.AD = coalesce(
        safe_get(d, paste0("count_cases_", assoc_chc_ad)),
        safe_get(d, paste0("count_cases_", assoc_cntr_ad))
      ),

      Freq.CHC = coalesce(
        safe_get(d, paste0("freq_controls_", assoc_chc_ad)),
        safe_get(d, paste0("freq_controls_", assoc_chc_cntr))
      ),
      Freq.CNTR = coalesce(
        safe_get(d, paste0("freq_cases_", assoc_chc_cntr)),
        safe_get(d, paste0("freq_controls_", assoc_cntr_ad))
      ),
      Freq.AD = coalesce(
        safe_get(d, paste0("freq_cases_", assoc_chc_ad)),
        safe_get(d, paste0("freq_cases_", assoc_cntr_ad))
      ),

      Beta.CHC_CNTR = safe_get(d, paste0("beta_", assoc_chc_cntr)),
      SE.CHC_CNTR = safe_get(d, paste0("standard_error_", assoc_chc_cntr)),
      P.CHC_CNTR = safe_get(d, paste0("p_value_", assoc_chc_cntr)),
      ORCI.CHC_CNTR = fmt_orci(
        safe_get(d, paste0("odds_ratio_", assoc_chc_cntr)),
        safe_get(d, paste0("lower_ci_", assoc_chc_cntr)),
        safe_get(d, paste0("upper_ci_", assoc_chc_cntr))
      ),

      Beta.CHC_AD = safe_get(d, paste0("beta_", assoc_chc_ad)),
      SE.CHC_AD = safe_get(d, paste0("standard_error_", assoc_chc_ad)),
      P.CHC_AD = safe_get(d, paste0("p_value_", assoc_chc_ad)),
      ORCI.CHC_AD = fmt_orci(
        safe_get(d, paste0("odds_ratio_", assoc_chc_ad)),
        safe_get(d, paste0("lower_ci_", assoc_chc_ad)),
        safe_get(d, paste0("upper_ci_", assoc_chc_ad))
      ),

      Beta.CNTR_AD = safe_get(d, paste0("beta_", assoc_cntr_ad)),
      SE.CNTR_AD = safe_get(d, paste0("standard_error_", assoc_cntr_ad)),
      P.CNTR_AD = safe_get(d, paste0("p_value_", assoc_cntr_ad)),
      ORCI.CNTR_AD = fmt_orci(
        safe_get(d, paste0("odds_ratio_", assoc_cntr_ad)),
        safe_get(d, paste0("lower_ci_", assoc_cntr_ad)),
        safe_get(d, paste0("upper_ci_", assoc_cntr_ad))
      )
    )
}

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

# Regression results for HLA-II haplotype genetic features (lead SNP, major
# MHC-II structure, one-field and two-field HLA alleles) across CHC, control
# and AD groups
regression_results <- fread("data_prep/study_cohort_genetics/output/hla_alleles_snps_regression.csv")

# Recode lead-SNP dosage columns to their rs-number for display; DR-broad
# major MHC-II structures are labeled DR1/DR2/DR4 directly in the input.
lead_snp_rs <- c(
  "chr6:32592593:G:T_T" = "rs35472547", # Hap-B lead SNP
  "chr6:32447376:C:T_T" = "rs9469112",  # Hap-R lead SNP
  "chr6:32609869:A:G_A" = "rs4335021"   # Hap-Y lead SNP
)
regression_results$allele <- ifelse(
  regression_results$allele %in% names(lead_snp_rs),
  unname(lead_snp_rs[regression_results$allele]),
  regression_results$allele
)

## -----------------------------------------------------------------------------
## 2. Define the genetic features of each haplotype, across resolution levels
## -----------------------------------------------------------------------------

# Each haplotype is represented at four resolution levels: lead SNP (rs),
# major MHC-II structure (one-field DR code and its underlying one-field
# DRB1/DQA1/DQB1), and two-field HLA-DRB1/DQA1/DQB1 alleles
hap_b_features <- c("rs35472547", "DR4", "DRB1*04", "DRB1*04:01", "DQA1*03", "DQA1*03:01", "DQB1*03", "DQB1*03:02")
hap_r_features <- c("rs9469112", "DR1", "DRB1*01", "DRB1*01:01", "DQA1*01", "DQA1*01:01", "DQB1*05", "DQB1*05:01")
hap_y_features <- c("rs4335021", "DR2", "DRB1*15", "DRB1*15:01", "DQA1*01", "DQA1*01:02", "DQB1*06", "DQB1*06:02")

## -----------------------------------------------------------------------------
## 3. Build Supplementary Table 1
## -----------------------------------------------------------------------------

supptable1 <- bind_rows(
  build_haplotype_table(regression_results, "Hap-B", hap_b_features),
  build_haplotype_table(regression_results, "Hap-R", hap_r_features),
  build_haplotype_table(regression_results, "Hap-Y", hap_y_features)
)

fwrite(supptable1, "tables/supptable1.tsv", sep = "\t")

## -----------------------------------------------------------------------------
## 4. Round Supplementary Table 1 for readability, and derive Table 2
## -----------------------------------------------------------------------------

supptable1_rounded <- supptable1 %>%
  mutate(
    across(c(N.CHC, N.CNTR, N.AD, Count.CHC, Count.CNTR, Count.AD), ~ round(.x, 0)),
    across(c(Freq.CHC, Freq.CNTR, Freq.AD), ~ round(.x, 4)),
    across(c(Beta.CHC_CNTR, SE.CHC_CNTR, Beta.CHC_AD, SE.CHC_AD, Beta.CNTR_AD, SE.CNTR_AD), ~ round(.x, 4))
  )

# Table 2: features reaching nominal significance (P<0.05) in at least one of
# the three pairwise comparisons
table2 <- supptable1_rounded %>% filter(P.CHC_CNTR < 0.05 | P.CHC_AD < 0.05 | P.CNTR_AD < 0.05)

supptable1_rounded <- supptable1_rounded %>%
  mutate(
    P.CHC_CNTR = formatC(P.CHC_CNTR, format = "e", digits = 2),
    P.CHC_AD = formatC(P.CHC_AD, format = "e", digits = 2),
    P.CNTR_AD = formatC(P.CNTR_AD, format = "e", digits = 2)
  )

fwrite(supptable1_rounded, "tables/supptable1_rounded.tsv", sep = "\t")

table2 <- table2 %>%
  mutate(
    P.CHC_CNTR = formatC(P.CHC_CNTR, format = "e", digits = 2),
    P.CHC_AD = formatC(P.CHC_AD, format = "e", digits = 2),
    P.CNTR_AD = formatC(P.CNTR_AD, format = "e", digits = 2)
  )
fwrite(table2, "tables/table2.tsv", sep = "\t")
