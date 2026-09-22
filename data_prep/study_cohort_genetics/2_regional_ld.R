## =============================================================================
## Regional LD (Figure 1): r2 of every SNP in the HLA/MHC region against the
## two EADB-GWAS-2026 top SNPs
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   EADB_GWAS_SUMSTATS_FULL (raw_input_data/data_paths.R)
##   data_prep/study_cohort_genetics/output/study_cohort_samples.txt - PARTICIPANT IDs, never
##     committed (see .gitignore) - from 1_clumping_cojo.R
##   data_prep/study_cohort_genetics/output/snps/chr6_mhc.dose.unscrambled.{pgen,pvar,psam} - from 1_clumping_cojo.R
## Outputs:
##   data_prep/study_cohort_genetics/output/regional_ld/regional_ld_pairs.ld - r2 of every SNP in the region
##     against the two EADB-GWAS-2026 top SNPs (one row per SNP pair), read directly by figure1_figure2.R
## =============================================================================

library(data.table)

# Run from the repository root, e.g. `Rscript data_prep/study_cohort_genetics/2_regional_ld.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

OUTDIR <- "data_prep/study_cohort_genetics/output/regional_ld"
STUDY_COHORT_SAMPLES <- "data_prep/study_cohort_genetics/output/study_cohort_samples.txt"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. Restrict the GWAS summary statistics to the HLA/MHC region
## -----------------------------------------------------------------------------

CHR <- 6

PROTECTIVE_TOP_SNP <- "6:32592593:G:T" # EADB-GWAS-2026 top SNP, rs35472547
RISK_TOP_SNP       <- "6:32447376:C:T" # EADB-GWAS-2026 top SNP, rs9469112

fullstats <- fread(EADB_GWAS_SUMSTATS_FULL, h = T, stringsAsFactors = F)

# Region spans the HLA class II genes, bounded by the first SNP at or after
# 32 Mb and the last SNP at or before 34 Mb on chromosome 6 - must match
# figure1_figure2.R's plot_start/plot_end exactly, since this script's LD
# matrix is plotted over that same region in Figure 1
REGION_START <- min(fullstats$POS[fullstats[["#CHROM"]] == CHR & fullstats$POS >= 32000000])
REGION_END   <- max(fullstats$POS[fullstats[["#CHROM"]] == CHR & fullstats$POS >= 33000000 & fullstats$POS <= 34000000])

fullstats_mhc <- fullstats[get("#CHROM") == CHR & get("POS") >= REGION_START & get("POS") <= REGION_END]

# The region can contain multiple rows per SNP ID (differing ref/alt allele
# coding); keep one row per ID
message("SNP IDs with more than one row in the region: ", sum(duplicated(fullstats_mhc$ID)))
fullstats_mhc <- unique(fullstats_mhc)
fullstats_mhc <- fullstats_mhc[!duplicated(ID)]

write.table(fullstats_mhc$ID, file = file.path(OUTDIR, "fullstats_gwas_snps_for_ld_matrix.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

## -----------------------------------------------------------------------------
## 2. Restrict cohort genotypes to the same samples and region, and compute LD
##    against the two top SNPs
## -----------------------------------------------------------------------------

system(sprintf(
  '%s --make-bed --bfile data_prep/study_cohort_genetics/output/snps/chr6_mhc.dose.unscrambled --keep %s --out %s/samples_kept',
  PLINK2_PATH, STUDY_COHORT_SAMPLES, OUTDIR))

plink_snps <- fread(file.path(OUTDIR, "samples_kept.bim"), header = FALSE)
setnames(plink_snps, old = c("V1", "V2"), new = c("CHR", "ID"), skip_absent = TRUE)
message("SNPs in PLINK file but not in GWAS summary stats: ", length(setdiff(plink_snps$ID, fullstats_mhc$ID)))
message("SNPs in GWAS summary stats but not in PLINK file: ", length(setdiff(fullstats_mhc$ID, plink_snps$ID)))

system(sprintf(
  '%s --make-bed --bfile %s/samples_kept --out %s/regional_ld_pairs --extract %s/fullstats_gwas_snps_for_ld_matrix.txt',
  PLINK2_PATH, OUTDIR, OUTDIR, OUTDIR))

writeLines(c(PROTECTIVE_TOP_SNP, RISK_TOP_SNP), file.path(OUTDIR, "ld_target_snps.txt"))
system(sprintf(
  '%s --bfile %s/regional_ld_pairs --r2 --ld-snp-list %s/ld_target_snps.txt --ld-window-r2 0 --ld-window-kb 5000 --ld-window 999999 --out %s/regional_ld_pairs',
  PLINK1_PATH, OUTDIR, OUTDIR, OUTDIR))

# Keep only the final pairwise LD file; drop the PLINK working files
file.remove(Sys.glob(file.path(OUTDIR, "samples_kept.*")))
file.remove(Sys.glob(file.path(OUTDIR, paste0("regional_ld_pairs.", c("bed", "bim", "fam", "log", "nosex")))))
file.remove(file.path(OUTDIR, c("fullstats_gwas_snps_for_ld_matrix.txt", "ld_target_snps.txt")))
