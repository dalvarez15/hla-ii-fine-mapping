## =============================================================================
## Combine HLA-II allele dosages with the literature and COJO haplotype SNP
## dosages, minor-allele oriented, plus DR-broad and one-field allele columns
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   data_prep/study_cohort_genetics/output/snps/chr6_mhc.raw - PLINK dosage export for
##     the HLA/MHC-region SNPs (from 1_clumping_cojo.R; ~5.7GB, not mirrored in raw_input_data/)
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - from 3_cojo_annotate.R
##   STUDY_COHORT_PHENOTYPES (raw_input_data/data_paths.R) - participant phenotypes, diagnosis,
##     sex, age and genetic PCs; no generating script in this repository
##   data_prep/study_cohort_genetics/output/hla_imputation/hla_allele_dosages_study_cohort.csv - from 4_read_hla_imputation.R
##   STUDY_COHORT_HLA_IMPUTATION_DIR (raw_input_data/data_paths.R) - re-read directly (not via
##     4_read_hla_imputation.R's output) to report imputation QC stats for this script's final cohort
## Outputs:
##   data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv - combined
##     SNP + two-field + one-field + DR-broad allele dosages; consumed by 6_haplotype_ld.R,
##     7_allele_snp_regression.R and ../100plus_study_cohort/1_combine_snps_alleles_neuropathology_100plus.R.
##     6_haplotype_ld.R selects specific columns for figure3.R's two-field-only LD matrix - see
##     that script for exactly which columns and why.
## =============================================================================

library(data.table)

# Run from the repository root, e.g. `Rscript data_prep/study_cohort_genetics/5_get_snps_alleles_dosages.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

dir.create("data_prep/study_cohort_genetics/output", recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Compares the allele frequency computed from `dosage_cols` in `data` against
# `expected_freq`, per SNP; used to confirm that allele-flip operations below
# produced the expected (minor-allele) orientation.
check_allele_frequencies <- function(data, snps, dosage_cols, expected_freq, tolerance = 0.01) {
    calculated <- sapply(dosage_cols, function(col) {
        if (is.na(col) || !(col %in% colnames(data))) return(NA_real_)
        vals <- data[[col]]
        sum(vals, na.rm = TRUE) / (sum(!is.na(vals)) * 2)
    })
    mismatched <- which(abs(calculated - expected_freq) > tolerance & !is.na(calculated))
    if (length(mismatched) > 0) {
        stop("Allele frequency mismatch after flipping for: ", paste(snps[mismatched], collapse = ", "))
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## 1. Select literature and COJO haplotype SNPs from the cohort's SNP dosages
## -----------------------------------------------------------------------------

SNPS_RAW_PATH <- "data_prep/study_cohort_genetics/output/snps/chr6_mhc.raw"

# The full file is ~5.7GB with 150k+ SNP columns, but only a handful of SNPs
# are actually needed below; reading the header alone (cheap) first lets us
# fread() only the needed columns afterwards, rather than loading everything
# into memory just to subset it away.
snps_header_chr <- paste0("chr", colnames(fread(SNPS_RAW_PATH, nrows = 0)))

# Previously reported AD/longevity SNPs (see figure1_figure2.R's previous_studies_lookup)
rs6605556   <- grep("chr6:32615322:A:G_A", snps_header_chr, value = TRUE) # Bellenguez 2022
rs9271192   <- grep("chr6:32610753:C:A_C", snps_header_chr, value = TRUE) # Lambert 2013
rs601945    <- grep("chr6:32605638:A:G_A", snps_header_chr, value = TRUE) # Le Guen shared AD-PD 2024
rs35472547  <- grep("chr6:32592593:G:T_G", snps_header_chr, value = TRUE) # Le Guen AD 2024 / EADB-GWAS-2026 protective lead SNP
rs34831921  <- grep("chr6:32622991:C:A_C", snps_header_chr, value = TRUE) # Joshi longevity 2017
rs9275152   <- grep("chr6:32684419:T:C_T", snps_header_chr, value = TRUE) # Timmers parental longevity 2019
hla_dra1    <- grep("chr6:32411770", snps_header_chr, value = TRUE)       # APOE-stratified GWAS, HLA-DRA SNP 1
hla_dra2    <- grep("chr6:32464090", snps_header_chr, value = TRUE)       # APOE-stratified GWAS, HLA-DRA SNP 2
rs_kunkle   <- grep("chr6:32607629", snps_header_chr, value = TRUE)       # Kunkle 2019

# COJO haplotype lead SNPs (minor-allele oriented)
cojo_ann <- read.csv2("data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv", header = TRUE, check.names = FALSE)
matches_snps <- unname(sapply(cojo_ann$SNP, function(snp) grep(snp, snps_header_chr, value = TRUE)))

snps_to_select <- unique(c(rs6605556, rs9271192, rs601945, rs35472547, rs34831921, rs9275152, matches_snps, hla_dra1, hla_dra2, rs_kunkle))
col_indices <- c(2, which(snps_header_chr %in% snps_to_select)) # column 2 = IID

selected_snps <- fread(SNPS_RAW_PATH, select = col_indices, header = TRUE, check.names = FALSE, data.table = FALSE)
colnames(selected_snps) <- paste0("chr", colnames(selected_snps))
colnames(selected_snps)[1] <- "ID_GWAS"

## -----------------------------------------------------------------------------
## 2. Merge phenotypes, keep participants passing QC in the study groups
## -----------------------------------------------------------------------------

phenotypes <- read.delim(STUDY_COHORT_PHENOTYPES)
selected_snps <- merge(selected_snps, phenotypes, by = "ID_GWAS", all.x = TRUE, all.y = FALSE)
selected_snps <- selected_snps[!is.na(selected_snps$PC1), ] # passed genetic-PC QC

our_phenos <- c("Centenarian", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD", "Control_MS", "Probable_AD", "AD_path", "Possible AD")
selected_snps_phenos <- subset(selected_snps, Diagnosis %in% our_phenos)

## -----------------------------------------------------------------------------
## 3. Re-orient COJO SNPs to the minor allele (matching cojo_ann), and flip
##    the literature SNPs to their minor allele too
## -----------------------------------------------------------------------------

for (i in seq_along(cojo_ann$SNP)) {
    snp <- cojo_ann$SNP[i]
    refA <- cojo_ann$refA[i]
    altA <- cojo_ann$altA[i]
    colname <- grep(snp, colnames(selected_snps_phenos), value = TRUE, ignore.case = FALSE)
    if (length(colname) == 1) {
        allele <- sub(".*_([A-Z])$", "\\1", colname)
        if (allele == refA) {
            next # already correct
        } else if (allele == altA) {
            selected_snps_phenos[[colname]] <- 2 - selected_snps_phenos[[colname]]
            new_colname <- sub(paste0("_", altA, "$"), paste0("_", refA), colname)
            setnames(selected_snps_phenos, colname, new_colname)
        } else {
            warning(sprintf("SNP %s column %s allele %s does not match refA (%s) or altA (%s)", snp, colname, allele, refA, altA))
        }
    } else if (length(colname) > 1) {
        warning(sprintf("Multiple columns found for SNP %s: %s", snp, paste(colname, collapse = ", ")))
    } else {
        warning(sprintf("No column found for SNP %s", snp))
    }
}
check_allele_frequencies(
    selected_snps_phenos, cojo_ann$SNP,
    sapply(cojo_ann$SNP, function(snp) { cn <- grep(snp, colnames(selected_snps_phenos), value = TRUE); if (length(cn) == 1) cn else NA }),
    cojo_ann$freq_geno
)

# Flip the literature SNPs (not part of the COJO table) to their minor allele
snp_cols_to_check <- c(rs6605556, rs9271192, rs601945, rs35472547, rs34831921, rs9275152, hla_dra1, hla_dra2, rs_kunkle)
for (col in snp_cols_to_check) {
    if (!(col %in% colnames(selected_snps_phenos))) next
    vals <- selected_snps_phenos[[col]]
    freq <- sum(vals, na.rm = TRUE) / (sum(!is.na(vals)) * 2)
    if (freq > 0.5) {
        selected_snps_phenos[[col]] <- 2 - vals
        allele_match <- regmatches(col, regexec("_(.)$", col))[[1]]
        allele <- if (length(allele_match) > 1) allele_match[2] else NA
        alleles <- unlist(strsplit(sub(".*:([ACGT]):([ACGT])_.*", "\\1:\\2", col), ":"))
        new_col <- col
        if (!is.na(allele) && length(alleles) == 2) {
            if (allele == alleles[1]) {
                new_col <- sub(paste0("_", alleles[1], "$"), paste0("_", alleles[2]), col)
            } else if (allele == alleles[2]) {
                new_col <- sub(paste0("_", alleles[2], "$"), paste0("_", alleles[1]), col)
            }
        }
        setnames(selected_snps_phenos, col, new_col)
    }
}

## -----------------------------------------------------------------------------
## 4. Merge with HLA allele dosages, derive DR-broad ancestral structures
## -----------------------------------------------------------------------------

alleles <- fread('data_prep/study_cohort_genetics/output/hla_imputation/hla_allele_dosages_study_cohort.csv', h = T, sep = ',', stringsAsFactors = F)
alleles_snps <- merge(alleles, selected_snps_phenos, by = "ID_GWAS")

# HLA imputation QC, restricted to this final cohort (cited in the manuscript).
# Posterior probability threshold as recommended by HIBAG's authors (Zheng et
# al. 2014, Pharmacogenomics J.) - see 4_read_hla_imputation.R, which applied it.
HIBAG_ACCEPT_PROB <- 0.5
imputation_raw <- do.call(rbind, lapply(c("DQA1", "DQB1", "DRB1"), function(locus) {
    read.delim(paste0(STUDY_COHORT_HLA_IMPUTATION_DIR, "/result_", locus, ".txt"))
}))
imputation_final_cohort <- imputation_raw[imputation_raw$sample.id %in% alleles_snps$ID_GWAS, ]
n_accepted <- sum(imputation_final_cohort$prob > HIBAG_ACCEPT_PROB)
n_total <- nrow(imputation_final_cohort)
message(sprintf(
    "HLA imputation QC: %d of %d (%.1f%%) imputed HLA genotypes accepted (posterior probability > %.1f) across HLA-DRB1, HLA-DQA1 and HLA-DQB1 in %d individuals",
    n_accepted, n_total, 100 * n_accepted / n_total, HIBAG_ACCEPT_PROB, length(unique(imputation_final_cohort$sample.id))
))

# Ancestral major MHC-II structures, from two-field HLA-DRB1 alleles
alleles_snps[, "DR1" := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns("DRB1\\*01|DRB1\\*10")]
alleles_snps[, "DR8" := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns("DRB1\\*08")]             # derivative of DR3
alleles_snps[, "DR2" := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns("DRB1\\*15|DRB1\\*16")]
alleles_snps[, "DR3" := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns("DRB1\\*03|DRB1\\*11|DRB1\\*12|DRB1\\*13|DRB1\\*14")]
alleles_snps[, "DR4" := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns("DRB1\\*04|DRB1\\*07|DRB1\\*09")]
alleles_snps[, "DRB1*04" := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns("DRB1\\*04")] # one-field DRB1*04, kept separately from the DR4 structure above

haplotype_cols <- c("DR1", "DR8", "DR2", "DR3", "DR4", "DRB1*04")
alleles_snps[, (haplotype_cols) := lapply(.SD, as.integer), .SDcols = haplotype_cols]
stopifnot(all(sapply(haplotype_cols, function(col) all(alleles_snps[[col]] %in% c(0L, 1L, 2L, NA)))))

# Round all SNP and two-field allele dosage columns to the nearest integer,
# before deriving one-field columns below so they sum already-rounded values
# (matches this script's behavior prior to merging into a single output file)
dosage_cols <- grep("^chr", colnames(alleles_snps), value = TRUE)
alleles_snps[, (dosage_cols) := lapply(.SD, function(x) as.integer(round(x))), .SDcols = dosage_cols]

## -----------------------------------------------------------------------------
## 5. Add one-field allele columns (e.g. "DRB1*01"), for 7_allele_snp_regression.R
## -----------------------------------------------------------------------------

for (f in c("01", "10", "08", "15", "16", "03", "11", "12", "13", "14", "07", "09")) {
    new_col <- paste0("DRB1*", f)
    alleles_snps[, (new_col) := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns(paste0("DRB1\\*", f))]
}

for (g in c("DQB1", "DQA1")) {
    g_cols <- grep(paste0("^", g, "\\*"), colnames(alleles_snps), value = TRUE)
    g_first <- unique(sub(paste0("^", g, "\\*([0-9]+).*"), "\\1", g_cols))
    for (f in g_first) {
        new_col <- paste0(g, "*", f)
        patt <- paste0(g, "\\*", f, "(:|$)") # ensure distinct matching (e.g. "03" doesn't also match "03:19")
        alleles_snps[, (new_col) := ifelse(rowSums(!is.na(.SD)) == 0, NA, rowSums(.SD, na.rm = TRUE)), .SDcols = patterns(patt)]
        alleles_snps[, (new_col) := as.integer(get(new_col))]
    }
}

write.csv2(alleles_snps, file = "data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv", row.names = FALSE)
