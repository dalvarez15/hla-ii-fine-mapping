## =============================================================================
## Supplementary Table 2 & Table 3: Hap-B/Hap-R/Hap-Y vs. microglial markers
## and neuropathology in the 100-plus Study and NBB replication cohort
##
## Formats the full association statistics into Supplementary Table 2, and
## FDR-significant subset into Table 3. The regressions themselves are fit on
## private individual-level data by data_prep/nbb_replication_cohort/3_microglia_neuropathology_regression.R;
## this script only reads that script's already-aggregate output, so unlike
## it, this one needs no private data. Figure 5 is likewise built from
## Supplementary Table 2's output by the separate figure5.R.
## =============================================================================
##
## Inputs:
##   data_prep/nbb_replication_cohort/output/microglia_neuropathology_regression.csv - full
##     association statistics (100-plus + NBB); from
##     data_prep/nbb_replication_cohort/3_microglia_neuropathology_regression.R
## Outputs:
##   tables/supptable2.tsv - read by figure5.R
##   tables/table3.tsv
## =============================================================================

library(data.table)
library(dplyr)

# Run from the repository root, e.g. `Rscript table3_supptable2.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## 1. Load the full association statistics and write Supplementary Table 2
## -----------------------------------------------------------------------------

supptable2 <- fread("data_prep/nbb_replication_cohort/output/microglia_neuropathology_regression.csv")

fwrite(
  supptable2,
  file = "tables/supptable2.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

## -----------------------------------------------------------------------------
## 2. Table 3: FDR-significant associations only, with rounded values
## -----------------------------------------------------------------------------

table3 <- supptable2 %>%
  filter(`FDR P` < 0.05) %>%
  mutate(
    P = formatC(P, format = "e", digits = 2),
    `FDR P` = formatC(`FDR P`, format = "e", digits = 2),
    Beta = round(Beta, 4),
    SE = round(SE, 4),
    `Std. Beta` = round(`Std. Beta`, 4)
  )

fwrite(table3, "tables/table3.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
