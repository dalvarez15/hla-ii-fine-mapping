## =============================================================================
## Figure 4 & Supplementary Figure 1: allele frequencies of Hap-B, Hap-R and
## Hap-Y (lead SNP, DR-broad code, two-field HLA-DRB1/DQA1/DQB1 alleles)
## across CHC, control (CNTR) and AD groups, with pairwise odds ratios and
## significance (logistic regression, corrected for sex and genetic PCs 1-5).
## Figure 4 shows the HLA-DRB1 locus only; Supplementary Figure 1 shows the
## full haplotype resolution.
## =============================================================================
##
## Inputs:
##   data_prep/study_cohort_genetics/output/hla_alleles_snps_regression.csv - regression
##     results for HLA-II haplotype genetic features: allele frequencies, sample sizes, and
##     CHC/control/AD association statistics; from 7_allele_snp_regression.R
## Outputs:
##   figures/figure4_R.pdf
##   figures/figuresupp1_R.pdf
## =============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Run from the repository root, e.g. `Rscript figure4_suppfigure1.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Display labels for the lead SNPs, by rs number, used throughout the plots
# below in place of the internal dosage-column names.
allele_display_labels <- c(
  "chr6:32592593:G:T_T" = "rs35472547", # Hap-B lead SNP
  "chr6:32447376:C:T_T" = "rs9469112",  # Hap-R lead SNP
  "chr6:32609869:A:G_A" = "rs4335021"   # Hap-Y lead SNP
)

recode_hla_labels <- function(x) {
  ifelse(x %in% names(allele_display_labels), unname(allele_display_labels[x]), x)
}

# Builds a bar plot of allele frequency (%) in CHCs, controls (CNTR) and AD
# patients for one haplotype's genetic features, with pairwise odds ratios
# and significance brackets (CHC vs. CNTR, CHC vs. AD, CNTR vs. AD).
#
# Input:
#   assocs           - the three association labels (in `regression`) for the
#                        CHC-vs-AD, CHC-vs-control and control-vs-AD models
#   features         - genetic features (alleles/SNPs) of the haplotype to plot,
#                        in display order
#   small_x_text     - if TRUE, shrink the x-axis facet labels (used for the
#                        full-resolution Supplementary Figure 1 panels, which
#                        have longer two-field allele names)
# Output: a ggplot object, one facet per feature
plot_allele_frequency_comparison <- function(assocs, features, small_x_text = FALSE) {
  regs <- regression %>%
    filter(allele %in% features,
           association %in% assocs) %>%
    select(allele, association, freq_cases, freq_controls, p_value, num_cases, num_controls, odds_ratio, lower_ci, upper_ci)
  if (nrow(regs) == 0) {
    message("No rows for associations: ", paste(assocs, collapse = ", "), " -> skipping")
    return(invisible(NULL))
  }

  regs$allele <- recode_hla_labels(regs$allele)
  features <- recode_hla_labels(features)

  wide <- regs %>%
    pivot_wider(
      names_from = association,
      values_from = c(freq_cases, freq_controls, p_value, num_cases, num_controls, odds_ratio, lower_ci, upper_ci),
      names_sep = "_"
    )

  a1 <- assocs[1] # CHC vs. AD
  a2 <- assocs[2] # CHC vs. CNTR
  a3 <- assocs[3] # CNTR vs. AD

  base_df <- data.frame(allele = features, stringsAsFactors = FALSE)
  joined <- base_df %>% left_join(wide, by = "allele")

  # CHC/CNTR/AD frequencies and pairwise stats are each reported under two of
  # the three association models (as "cases" or "controls" depending on which
  # comparison); coalesce() picks whichever model reports the value.
  freqs <- data.frame(
    allele = joined$allele,
    CHCs = coalesce(joined[[paste0("freq_controls_", a2)]], joined[[paste0("freq_controls_", a1)]]),
    CNTR = coalesce(joined[[paste0("freq_cases_", a2)]],   joined[[paste0("freq_controls_", a3)]]),
    AD   = coalesce(joined[[paste0("freq_cases_", a1)]],   joined[[paste0("freq_cases_", a3)]]),
    p_cent_ctrl  = joined[[paste0("p_value_", a2)]],
    p_cent_ad    = joined[[paste0("p_value_", a1)]],
    p_ctrl_ad    = joined[[paste0("p_value_", a3)]],
    or_cent_ctrl  = joined[[paste0("odds_ratio_", a2)]],
    or_cent_ad    = joined[[paste0("odds_ratio_", a1)]],
    or_ctrl_ad    = joined[[paste0("odds_ratio_", a3)]],
    lci_cent_ctrl = joined[[paste0("lower_ci_", a2)]],
    lci_cent_ad   = joined[[paste0("lower_ci_", a1)]],
    lci_ctrl_ad   = joined[[paste0("lower_ci_", a3)]],
    uci_cent_ctrl = joined[[paste0("upper_ci_", a2)]],
    uci_cent_ad   = joined[[paste0("upper_ci_", a1)]],
    uci_ctrl_ad   = joined[[paste0("upper_ci_", a3)]],
    n_cent      = coalesce(joined[[paste0("num_controls_", a2)]], joined[[paste0("num_controls_", a1)]]),
    n_controls  = coalesce(joined[[paste0("num_cases_", a2)]],    joined[[paste0("num_controls_", a3)]]),
    n_ad        = coalesce(joined[[paste0("num_cases_", a1)]],    joined[[paste0("num_cases_", a3)]]),
    stringsAsFactors = FALSE
  )

  freqs$allele <- factor(freqs$allele, levels = features)

  pick_first <- function(x) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0) return(NA_integer_)
    x2[1]
  }
  n_cent_val <- pick_first(freqs$n_cent)
  n_ctrl_val <- pick_first(freqs$n_controls)
  n_ad_val   <- pick_first(freqs$n_ad)

  x_labels <- c(
    "CHCs" = paste0("CHCs\nn=", ifelse(is.na(n_cent_val), "NA", n_cent_val)),
    "CNTR" = paste0("CNTR\nn=", ifelse(is.na(n_ctrl_val), "NA", n_ctrl_val)),
    "AD"   = paste0("AD\nn=",   ifelse(is.na(n_ad_val),   "NA", n_ad_val))
  )

  df_long <- freqs %>%
    pivot_longer(c(CHCs, CNTR, AD), names_to = "group", values_to = "freq") %>%
    mutate(group = factor(group, levels = c("CHCs", "CNTR", "AD")),
           freq_pct = 100 * freq)

  top <- df_long %>%
    group_by(allele) %>%
    summarise(ymax = max(freq_pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(ymax = ifelse(is.finite(ymax), ymax, 0))

  pv <- freqs %>%
    pivot_longer(cols = c(p_cent_ctrl, p_cent_ad, p_ctrl_ad),
                 names_to = "comp", values_to = "p") %>%
    mutate(
      group1 = case_when(
        comp == "p_cent_ctrl" ~ "CHCs",
        comp == "p_cent_ad"   ~ "CHCs",
        comp == "p_ctrl_ad"   ~ "CNTR"
      ),
      group2 = case_when(
        comp == "p_cent_ctrl" ~ "CNTR",
        comp == "p_cent_ad"   ~ "AD",
        comp == "p_ctrl_ad"   ~ "AD"
      ),
      or = case_when(
        comp == "p_cent_ctrl" ~ or_cent_ctrl,
        comp == "p_cent_ad"   ~ or_cent_ad,
        comp == "p_ctrl_ad"   ~ or_ctrl_ad
      ),
      lci = case_when(
        comp == "p_cent_ctrl" ~ lci_cent_ctrl,
        comp == "p_cent_ad"   ~ lci_cent_ad,
        comp == "p_ctrl_ad"   ~ lci_ctrl_ad
      ),
      uci = case_when(
        comp == "p_cent_ctrl" ~ uci_cent_ctrl,
        comp == "p_cent_ad"   ~ uci_cent_ad,
        comp == "p_ctrl_ad"   ~ uci_ctrl_ad
      ),
      label = case_when(
        is.na(p)  ~ "n.s.",
        p < 0.001 ~ paste0("***", sprintf("%.2f", or), " [", sprintf("%.2f", lci), "-", sprintf("%.2f", uci), "]"),
        p < 0.01  ~ paste0("**", sprintf("%.2f", or), " [", sprintf("%.2f", lci), "-", sprintf("%.2f", uci), "]"),
        p < 0.05  ~ paste0("*", sprintf("%.2f", or), " [", sprintf("%.2f", lci), "-", sprintf("%.2f", uci), "]"),
        TRUE      ~ "n.s."
      )
    ) %>%
    select(allele, group1, group2, p, label) %>%
    left_join(top, by = "allele") %>%
    group_by(allele) %>%
    mutate(
      step = pmax(3, 0.08 * (ymax + 1)),
      y.position = ymax + step * row_number()
    ) %>%
    ungroup()

  p <- ggplot(df_long, aes(group, freq_pct, fill = group)) +
    labs(x = NULL, y = "Allele frequency (%)") +
    geom_col(width = 0.6, color = "black") +
    facet_wrap(~ allele, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(
      "CHCs" = "#B8860B",
      "CNTR" = "#1F78B4",
      "AD"   = "#E41A1C"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    scale_x_discrete(labels = x_labels) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          strip.text = element_text(size = 10)) +
    geom_text(aes(label = ifelse(is.na(freq_pct), "", sprintf("%.1f", freq_pct))),
              vjust = -0.5, size = 3) +
    stat_pvalue_manual(
      pv,
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      label = "label",
      tip.length = 0.01,
      bracket.size = 0.3,
      size = 2.5,
      inherit.aes = FALSE
    )

  if (small_x_text) {
    p <- p + theme(axis.text.x = element_text(size = 7))
  }

  p
}

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

# Regression results for HLA-II haplotype genetic features across CHC,
# control (CNTR) and AD groups (CHC-vs-AD, CHC-vs-control, control-vs-AD)
regression <- fread("data_prep/study_cohort_genetics/output/hla_alleles_snps_regression.csv")

# The three pairwise comparisons plotted for every haplotype
comparison_associations <- c("cent_ad", "cent_ctrl", "ctrl_ad")

## -----------------------------------------------------------------------------
## 2. Identify each haplotype's genetic features, at DR-only and full
##    (lead SNP through HLA-DQB1) resolution
## -----------------------------------------------------------------------------

snp_targets <- c("Hap-B", "Hap-R", "Hap-Y")
snp_alleles <- regression %>%
  filter(snp_easy %in% snp_targets) %>%
  arrange(match(snp_easy, snp_targets)) %>% # keep order to display 
  pull(allele) %>%
  unique() %>%
  as.character()

#snp_alleles[1:2] <- snp_alleles[2:1] # swap so Hap-B's lead SNP comes first, matching haplotype display order

# Full haplotype resolution (lead SNP, DR-broad, one-field DRB1/DQA1/DQB1,
# two-field DRB1/DQA1/DQB1): Supplementary Figure 1
hap_b_features_full <- c(snp_alleles[1], "DR4", "DRB1*04", "DRB1*04:01", "DQA1*03", "DQA1*03:01", "DQB1*03", "DQB1*03:02")
hap_r_features_full <- c(snp_alleles[2], "DR1", "DRB1*01", "DRB1*01:01", "DQA1*01", "DQA1*01:01", "DQB1*05", "DQB1*05:01")
hap_y_features_full <- c(snp_alleles[3], "DR2", "DRB1*15", "DRB1*15:01", "DQA1*01", "DQA1*01:02", "DQB1*06", "DQB1*06:02")

# HLA-DRB1 locus only: Figure 4
hap_b_features_dr <- c(snp_alleles[1], "DR4", "DRB1*04", "DRB1*04:01")
hap_r_features_dr <- c(snp_alleles[2], "DR1", "DRB1*01", "DRB1*01:01")
hap_y_features_dr <- c(snp_alleles[3], "DR2", "DRB1*15", "DRB1*15:01")

## -----------------------------------------------------------------------------
## 3. Supplementary Figure 1: allele frequencies at full haplotype resolution
## -----------------------------------------------------------------------------

hap_b_plot_full <- plot_allele_frequency_comparison(comparison_associations, hap_b_features_full, small_x_text = TRUE)
hap_r_plot_full <- plot_allele_frequency_comparison(comparison_associations, hap_r_features_full, small_x_text = TRUE)
hap_y_plot_full <- plot_allele_frequency_comparison(comparison_associations, hap_y_features_full, small_x_text = TRUE)

suppfigure1 <- hap_b_plot_full / hap_r_plot_full / hap_y_plot_full
ggsave("figures/figuresupp1_R.pdf", suppfigure1, width = 14, height = 12)

## -----------------------------------------------------------------------------
## 4. Figure 4: allele frequencies at the HLA-DRB1 locus only
## -----------------------------------------------------------------------------

hap_b_plot_dr <- plot_allele_frequency_comparison(comparison_associations, hap_b_features_dr)
hap_r_plot_dr <- plot_allele_frequency_comparison(comparison_associations, hap_r_features_dr)
hap_y_plot_dr <- plot_allele_frequency_comparison(comparison_associations, hap_y_features_dr)

figure4 <- hap_b_plot_dr / hap_r_plot_dr / hap_y_plot_dr
ggsave("figures/figure4_R.pdf", figure4, width = 14, height = 12)
