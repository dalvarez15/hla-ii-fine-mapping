## =============================================================================
## LD-based clumping (r2=0.80) and COJO fine-mapping of the HLA/MHC region
## (chr6:28,510,120-33,480,575): tests whether clumps reaching P<=1e-5
## ("suggestive") and within 1 Mb of the region's EADB-GWAS-2026 top SNP
## (protective, rs35472547) remain independently associated with AD, via
## GCTA's Conditional Joint (COJO) analysis, then FDR-adjusts the joint
## p-values.
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   CHR6_GENOMEWIDE_DOSAGE, EADB_GWAS_SUMSTATS_FULL, STUDY_COHORT_PHENOTYPES (raw_input_data/data_paths.R)
## Outputs:
##   data_prep/study_cohort_genetics/output/study_cohort_samples.txt - PARTICIPANT IDs, never
##     committed (see .gitignore) - every participant with a diagnosis this analysis uses (same
##     filter as 5_get_snps_alleles_dosages.R), read directly by 2_regional_ld.R
##   data_prep/study_cohort_genetics/output/snps/chr6_mhc.dose.unscrambled.* - HLA/MHC-region dosage data
##   data_prep/study_cohort_genetics/output/snps/chr6_mhc.raw - per-sample dosage export, consumed by 5_get_snps_alleles_dosages.R
##   data_prep/study_cohort_genetics/output/clumping_r2_080.clumps - LD-based clumping output, read directly by figure1_figure2.R
##   data_prep/study_cohort_genetics/output/cojo/summary_cojo.txt - COJO results per lead SNP, consumed by 3_cojo_annotate.R
## =============================================================================

library(data.table)

# Run from the repository root, e.g. `Rscript data_prep/study_cohort_genetics/1_clumping_cojo.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

GENOTYPE_PATH <- "data_prep/study_cohort_genetics/output/snps/"
CLUMP_R2 <- 0.80
P_SUGGESTIVE <- 1e-5 # must match figure1_figure2.R's P_SUGGESTIVE

# EADB-GWAS-2026 top SNP for this locus (protective, rs35472547, P=2.61e-23);
# COJO tests clumps within 1 Mb of it for independent association
TOP_SNP_POS <- 32592593

dir.create("data_prep/study_cohort_genetics/output/snps", recursive = TRUE, showWarnings = FALSE)
dir.create("data_prep/study_cohort_genetics/output/cojo", recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. Define the study cohort: everyone who passed genetic-PC QC with a
##    diagnosis this analysis uses (same filter 5_get_snps_alleles_dosages.R
##    applies for its own cohort, verified to match it exactly)
## -----------------------------------------------------------------------------

STUDY_COHORT_SAMPLES <- "data_prep/study_cohort_genetics/output/study_cohort_samples.txt"
our_phenos <- c("Centenarian", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD", "Control_MS", "Probable_AD", "AD_path", "Possible AD")
phenotypes <- read.delim(STUDY_COHORT_PHENOTYPES)

# Cohort composition by recruitment source, as reported in the manuscript's
# cohort description - both before and after the genetic-PC QC filter this
# script applies below (`!is.na(PC1)`), which is what defines "our cohort"
# everywhere else in this pipeline. Counted by unique ID_GWAS, matching
# study_cohort_ids below: a small number of participants have more than one
# row in this file (e.g. genotyped on two chips) and must not be
# double-counted.
# "SCD" is how ADC's cognitively healthy controls are coded in Diagnosis;
# "Control_MS" is how Project Y's controls are coded (not multiple sclerosis
# controls, despite the label); NBB-recruited controls are coded
# "Control_path". Diagnosis alone doesn't distinguish NBB-recruited from
# ADC-recruited AD cases, since both can carry "AD_path" - that split uses
# Study (NBB's "hersenbank"/"Hersenbank" values).
control_sources <- c(
    LASA = "Control_LASA", `ADC (SCD)` = "SCD", `100-plus partners` = "Control_100plus",
    `Netherlands Twin Study` = "Control_other_twin", `Project Y (Control_MS)` = "Control_MS",
    `NBB (Control_path)` = "Control_path"
)
ad_diagnoses <- c("Probable_AD", "AD_path", "Possible AD")
ad_nbb <- phenotypes$Diagnosis == "AD_path" & phenotypes$Study %in% c("hersenbank", "Hersenbank")
not_duplicate_id <- !duplicated(phenotypes$ID_GWAS)

report_cohort_composition <- function(qc_pass, label) {
    qc_pass <- qc_pass & not_duplicate_id
    chc_n <- sum(phenotypes$Diagnosis == "Centenarian" & qc_pass)
    control_n <- sapply(control_sources, function(d) sum(phenotypes$Diagnosis == d & qc_pass))
    ad_n <- sum(phenotypes$Diagnosis %in% ad_diagnoses & qc_pass)
    ad_nbb_n <- sum(ad_nbb & qc_pass)
    message(sprintf("%s cohort composition: CHC=%d, Controls=%d (%s), AD=%d (NBB=%d, ADC/other=%d), Total=%d",
                     label, chc_n, sum(control_n), paste(sprintf("%s=%d", names(control_n), control_n), collapse = ", "),
                     ad_n, ad_nbb_n, ad_n - ad_nbb_n, chc_n + sum(control_n) + ad_n))
}
report_cohort_composition(rep(TRUE, nrow(phenotypes)), "Pre-QC")
report_cohort_composition(!is.na(phenotypes$PC1), "Post-QC (our cohort)")

study_cohort_ids <- unique(phenotypes$ID_GWAS[!is.na(phenotypes$PC1) & phenotypes$Diagnosis %in% our_phenos])
writeLines(as.character(study_cohort_ids), STUDY_COHORT_SAMPLES)
message("Study cohort: ", length(study_cohort_ids), " participants")

## -----------------------------------------------------------------------------
## 2. Extract the HLA/MHC region from the cohort's chr6 dosage data
## -----------------------------------------------------------------------------

system(paste0(
    PLINK2_PATH, " --pfile ", CHR6_GENOMEWIDE_DOSAGE, " ",
    "--chr 6 --from-bp 28510120 --to-bp 33480575 --make-pgen --out ", GENOTYPE_PATH, "chr6_mhc.dose.unscrambled"
))

# Per-sample dosage export (PLINK "A" format), used by 5_get_snps_alleles_dosages.R
system(paste0(
    PLINK2_PATH, " --pfile ", GENOTYPE_PATH, "chr6_mhc.dose.unscrambled --export A --out ", GENOTYPE_PATH, "chr6_mhc"
))

# PLINK 1 binary format (.bed/.bim/.fam), used by 2_regional_ld.R
system(paste0(
    PLINK2_PATH, " --pfile ", GENOTYPE_PATH, "chr6_mhc.dose.unscrambled --make-bed --out ", GENOTYPE_PATH, "chr6_mhc.dose.unscrambled"
))

## -----------------------------------------------------------------------------
## 3. LD-based clumping at r2=0.80, flag suggestive clumps
## -----------------------------------------------------------------------------

system(paste0(
    PLINK2_PATH, " --pfile ", GENOTYPE_PATH, "chr6_mhc.dose.unscrambled ",
    "--keep ", STUDY_COHORT_SAMPLES, " ",
    "--clump ", EADB_GWAS_SUMSTATS_FULL, " ",
    "--clump-r2 ", CLUMP_R2, " --clump-p1 0.01 --clump-p2 1 --out data_prep/study_cohort_genetics/output/clumping_r2_080"
))

# PLINK2 also writes a .clumps.missing_id listing every SNP in the genome-wide
# EADB_GWAS_SUMSTATS_FULL that isn't in this chr6-region-only genotype file -
# expected (this file is deliberately chr6-only) and unused downstream.
missing_id_file <- "data_prep/study_cohort_genetics/output/clumping_r2_080.clumps.missing_id"
if (file.exists(missing_id_file)) file.remove(missing_id_file)

clumps_df <- fread("data_prep/study_cohort_genetics/output/clumping_r2_080.clumps", h = T, stringsAsFactors = F)
clumps_df$suggestive <- ifelse(clumps_df$P <= P_SUGGESTIVE, 'yes', 'no')
message("Clumps reaching P<=", P_SUGGESTIVE, ": ", sum(clumps_df$suggestive == 'yes'), " of ", nrow(clumps_df))

## -----------------------------------------------------------------------------
## 4. Load chr6 HLA-II region specific summary statistics from EADB-GWAS-2026
## -----------------------------------------------------------------------------

fullstats <- fread(EADB_GWAS_SUMSTATS_FULL, h = T, stringsAsFactors = F)

CHR <- 6
tmp_sumstats <- fullstats[which(fullstats$"#CHROM" == CHR), ]
# Align effect allele with the genotyping data's ALT allele
tmp_sumstats$freq_checked <- ifelse((tmp_sumstats$beta_ALT * tmp_sumstats$beta) > 0, tmp_sumstats$effect_allele_frequency, 1 - tmp_sumstats$effect_allele_frequency)
tmp_sumstats$total_n <- tmp_sumstats$n_cases + tmp_sumstats$n_controls
tmp_sumstats <- tmp_sumstats[, c('ID', 'ALT', 'REF', 'freq_checked', 'beta_ALT', 'standard_error', 'P', 'total_n')]
colnames(tmp_sumstats) <- c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N')

## -----------------------------------------------------------------------------
## 5. COJO: test independence of suggestive clumps within 1 Mb of the top SNP
## -----------------------------------------------------------------------------

outdir <- "data_prep/study_cohort_genetics/output/cojo/"
tmp_clumps <- clumps_df[which(clumps_df$suggestive == 'yes'), ]
tmp_clumps_snps <- tmp_clumps[which(abs(tmp_clumps$POS - TOP_SNP_POS) <= 1000000), ]
tmp_fullstats_top <- tmp_sumstats[which(tmp_sumstats$SNP %in% tmp_clumps_snps$ID), ]

write.table(tmp_fullstats_top[!duplicated(tmp_fullstats_top$SNP), ], paste0(outdir, 'tmp_cojo_top.ma'), quote = F, row.names = F, sep = "\t")
write.table(tmp_fullstats_top$SNP, paste0(outdir, 'tmp_snps_top.txt'), quote = F, row.names = F, col.names = F)

# Make a small PLINK file for just these SNPs, for speed
cmd_plink <- paste0(PLINK2_PATH, " --pfile ", GENOTYPE_PATH, 'chr', CHR, '_mhc.dose.unscrambled --extract ', outdir, 'tmp_snps_top.txt --hard-call-threshold 0.4 --keep ', STUDY_COHORT_SAMPLES, ' --make-bed --out ', outdir, 'tmp_plink')
system(cmd_plink)

cmd_top <- paste0(GCTA_PATH, ' --bfile ', outdir, 'tmp_plink', ' --chr ', CHR, ' --maf 0.01 --cojo-file ', outdir, 'tmp_cojo_top.ma', ' --cojo-slct --out ', outdir, 'cojo_top --extract ', outdir, 'tmp_snps_top.txt --cojo-p 0.05 --threads 4')
system(cmd_top)
system(paste0('rm ', outdir, 'tmp*'))

## -----------------------------------------------------------------------------
## 6. FDR-correct the joint COJO p-values and write the final file
## -----------------------------------------------------------------------------

summary_data <- fread(paste0(outdir, 'cojo_top.jma.cojo'), h = TRUE, stringsAsFactors = FALSE)
summary_data <- summary_data[order(summary_data$pJ), ]
summary_data$pJ_FDR <- p.adjust(summary_data$pJ, method = 'fdr')
write.table(summary_data, paste0(outdir, 'summary_cojo.txt'), quote = F, row.names = F, sep = '\t')
