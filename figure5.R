## =============================================================================
## Figure 5: heatmap of standardized effect sizes (Hap-B/Hap-R/Hap-Y genetic
## features vs. microglial markers and neuropathology), 100-plus Study
## microglia panel, 100-plus Study neuropathology panel, and NBB replication
## cohort neuropathology panel, with FDR-significant associations labeled.
##
## Built entirely from Supplementary Table 2's already-aggregated association
## statistics (from table3_supptable2.R), so - like every other main script -
## it needs no private individual-level data and can be run from a clone of
## this repository as-is.
## =============================================================================
##
## Inputs:
##   tables/supptable2.tsv - full association statistics (Haplotype, Feature,
##     Outcome variable, Std. Beta, FDR P etc.); from table3_supptable2.R
## Outputs:
##   figures/figure5_R.pdf
## =============================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Run from the repository root, e.g. `Rscript figure5.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Heatmap of standardized effect sizes (feature x outcome), with FDR labels.
make_panel <- function(df, outcome_order, title, zlim) {
  df <- df %>%
    filter(test_var %in% outcome_order) %>%
    mutate(
      test_var = factor(test_var, levels = rev(outcome_order)),
      snp_var  = factor(snp_var, levels = haplotype_features)
    )

  ggplot(df, aes(x = snp_var, y = test_var, fill = beta_std)) +
    geom_tile(color = "grey80") +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-zlim, zlim), name = "Effect\n(std. beta)"
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

## -----------------------------------------------------------------------------
## 1. Panel groupings and haplotype feature order (matches
##    data_prep/nbb_replication_cohort/3_microglia_neuropathology_regression.R)
## -----------------------------------------------------------------------------

neuropathology_scores <- c("CERAD score", "Braak NFT stage (100+)", "Thal Aβ phase")
neuropathology_loads  <- c("Aβ load", "Aβ40 load", "Aβ42 load", "pTau217 load", "AT8 load", "GT-38 load")
microglia_loads       <- c("HLA class II load", "CD11c load", "CD68 load", "Iba1 load", "P2RY12 load")
neuropathology_all    <- c(neuropathology_scores, neuropathology_loads)

nbb_outcome_order <- c(
  "Braak amyloid-β stage",
  "Braak NFT stage (NBB)",
  "Braak α-synuclein stage",
  "Braak amyloid-β stage\ncorrected for:\nBraak NFT stage (NBB)",
  "Braak NFT stage (NBB)\ncorrected for:\nBraak amyloid-β stage"
)

haplotype_features <- c(
  "rs35472547", "DR4", "DRB1*04:01", "DQA1*03:01", "DQB1*03:02", # Hap-B (blue)
  "rs9469112",  "DR1", "DRB1*01:01", "DQA1*01:01", "DQB1*05:01", # Hap-R (red)
  "rs4335021",  "DR2", "DRB1*15:01", "DQA1*01:02", "DQB1*06:02"  # Hap-Y (yellow)
)

## -----------------------------------------------------------------------------
## 2. Rebuild the 100-plus and NBB association tables from Supplementary
##    Table 2, restoring the cross-adjustment labels collapsed for the TSV
## -----------------------------------------------------------------------------

supptable2 <- fread("tables/supptable2.tsv")
supptable2$`Outcome variable` <- gsub(" corrected for: ", "\ncorrected for:\n", supptable2$`Outcome variable`)

to_panel_df <- function(df) {
  data.frame(
    test_var = df$`Outcome variable`,
    snp_var  = factor(df$Feature, levels = haplotype_features),
    beta_std = df$`Std. Beta`,
    label    = ifelse(df$`FDR P` < 0.05, sprintf("fdr=%.2g", df$`FDR P`), "")
  )
}

results     <- to_panel_df(supptable2[supptable2$`Outcome variable` %in% c(neuropathology_all, microglia_loads), ])
nbb_results <- to_panel_df(supptable2[supptable2$`Outcome variable` %in% nbb_outcome_order, ])

## -----------------------------------------------------------------------------
## 3. Figure 5: heatmap panels (100-plus microglia, 100-plus neuropathology,
##    NBB neuropathology)
## -----------------------------------------------------------------------------

zlim <- max(abs(supptable2$`Std. Beta`), na.rm = TRUE)

plot_micro <- make_panel(results, microglia_loads, "Microglia 100+ (N.CHC=89, N.AD=7)", zlim)
plot_neuro <- make_panel(results, neuropathology_all, "Neuropathology 100+ (N.CHC=89, N.AD=7)", zlim)
plot_nbb <- make_panel(nbb_results, nbb_outcome_order, "Neuropathology NBB (N.CNTR=248, N.AD=513)", zlim)

cairo_pdf(
  "figures/figure5_R.pdf",
  width = 14,
  height = 17
)
spacer <- grid::rectGrob(gp = grid::gpar(col = NA, fill = NA)) # small vertical spacer between panels
grid.arrange(
  spacer, plot_micro,
  spacer, plot_neuro,
  spacer, plot_nbb,
  ncol = 1,
  heights = c(0.5, 4, 0.5, 7, 0.5, 4, 0.5)
)
dev.off()
