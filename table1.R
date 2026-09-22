## =============================================================================
## Table 1: the six COJO-independent HLA-II haplotypes (Hap-1, Hap-R, Hap-B,
## Hap-Y, Hap-5, Hap-6), each annotated with its major ancestral MHC-II
## structure, imputed two-field HLA-DRB1/DQA1/DQB1 alleles in LD (r2>0.2), and
## overlap with previously reported AD/longevity-associated SNPs (in LD, r2>0.2)
## =============================================================================
##
## Inputs:
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - COJO-identified
##     haplotype lead SNPs (minor-allele oriented), with GWAS and joint-model (COJO) association
##     statistics; from 3_cojo_annotate.R
## Outputs:
##   tables/table1.tsv
## =============================================================================

library(data.table)
library(dplyr)

# Run from the repository root, e.g. `Rscript table1.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

# COJO-independent haplotype lead SNPs (6 haplotypes in the HLA-II region),
# with GWAS and joint-model (COJO) association statistics
cojo_haplotypes <- read.csv2(
  "data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv"
)
cojo_haplotypes$ID <- sub("^chr", "", cojo_haplotypes$SNP)

# rs number of the three lead SNPs that overlap previously reported AD or
# longevity signals (Hap-R, Hap-B, Hap-Y; see step 2)
hap_rs <- c("Hap-B" = "rs35472547", "Hap-R" = "rs9469112", "Hap-Y" = "rs4335021")
cojo_haplotypes$rs <- unname(hap_rs[cojo_haplotypes$snp_easy])

# Lookup of SNPs previously reported in the literature as associated with AD
# or longevity, in LD (r2>0.2) with one of the six COJO haplotype lead SNPs.
# `ld` records which of our haplotypes (Hap-B/Hap-R/Hap-Y) each one tags.
previous_studies_lookup <- data.frame(
  ID = c(
    "NA", "NA",
    "6:32615322:A:G", "6:32610753:C:A", "6:32605638:A:G", "6:32607629:A:T",
    "6:32622991:C:A", "6:32684419:T:C",
    "6:32411770:C:T", "6:32464090:G:T"
  ),
  rs = c(
    "rs9469112", "rs35472547",
    "rs6605556", "rs9271192", "rs601945", "rs9271058",
    "rs34831921", "rs9275152",
    "rs17208902", "rs9268888"
  ),
  study = c(
    "EADB (2026)-1", "EADB (2026)-2",
    "Bellenguez (2022)", "Lambert (2013)", "Le Guen (2023)", "Kunkle (2019)",
    "Joshi (2017)", "Timmers (2019)",
    "Thomassen (2026)-1", "Thomassen (2026)-2"
  ),
  ld = c(
    "Hap-R", "Hap-B",
    "Hap-B", "Hap-Y", "Hap-B", "Hap-Y",
    "none", "Hap-B",
    "Hap-R", "Hap-Y"
  ),
  stringsAsFactors = FALSE
)
previous_studies_lookup$snp_study <- paste0(previous_studies_lookup$rs, ", ", previous_studies_lookup$study)

## -----------------------------------------------------------------------------
## 2. Build the haplotype table: GWAS and COJO association statistics
## -----------------------------------------------------------------------------

table1 <- cojo_haplotypes %>%
  mutate(Haplotype = snp_easy) %>%
  select(Haplotype, rs, SNP, Chr, bp, refA, altA, freq, b, se, p, bJ, bJ_se, pJ, pJ_FDR) %>%
  rename(
    `rs lead SNP` = rs,
    `ID lead SNP` = SNP,
    `Frequency GWAS` = freq,
    `Beta GWAS` = b,
    `SE GWAS` = se,
    `P GWAS` = p,
    `Beta COJO` = bJ,
    `SE COJO` = bJ_se,
    `P COJO` = pJ,
    `P COJO FDR` = pJ_FDR
  )

## -----------------------------------------------------------------------------
## 3. Annotate with previously reported signals and HLA-II alleles in LD
## -----------------------------------------------------------------------------

# Collapse the literature lookup to one "SNP, study" string per haplotype
previous_studies_by_hap <- previous_studies_lookup %>%
  filter(!is.na(ld), !is.na(snp_study), ld != "", snp_study != "") %>%
  group_by(ld) %>%
  summarise(
    `Previous associated SNPs in LD r2>0.2, Study` = paste(unique(snp_study), collapse = "; "),
    .groups = "drop"
  )

hap_dr_broad <- c("Hap-B" = "DR4", "Hap-R" = "DR1", "Hap-Y" = "DR2")
table1 <- table1 %>%
  left_join(previous_studies_by_hap, by = c("Haplotype" = "ld")) %>%
  mutate(
    `Major MHC-II structures in LD r2>0.2` = unname(hap_dr_broad[Haplotype]),
    `HLA-II alleles in LD r2>0.2` = case_when(
      Haplotype == "Hap-B" ~ "DRB1*04:04, DQA1*03:01, DQB1*03:02",
      Haplotype == "Hap-R" ~ "DRB1*01:01, DQA1*01:01, DQB1*05:01",
      Haplotype == "Hap-Y" ~ "DRB1*15:01, DQA1*01:02, DQB1*06:02",
      TRUE ~ NA_character_
    )
  )

## -----------------------------------------------------------------------------
## 4. Format and save Table 1
## -----------------------------------------------------------------------------

table1 <- table1 %>%
  mutate(
    `P GWAS` = formatC(`P GWAS`, format = "e", digits = 2),
    `P COJO` = formatC(`P COJO`, format = "e", digits = 2),
    `P COJO FDR` = formatC(`P COJO FDR`, format = "e", digits = 2)
  )

fwrite(
  table1,
  "tables/table1.tsv",
  sep = "\t"
)
