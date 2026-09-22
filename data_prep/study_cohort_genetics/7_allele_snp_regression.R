## =============================================================================
## HLA-II allele/SNP association regression (dosage ~ allele + PC1-5, additive
## genetic model, logistic) for the three group comparisons reported in the
## paper: control-vs-AD, CHC-vs-control, CHC-vs-AD.
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv - per-participant
##     HLA-II allele/SNP dosages (one-field and two-field), diagnosis, sex and age (from
##     5_get_snps_alleles_dosages.R)
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - COJO haplotype
##     lead SNPs with easy-to-read labels (from 3_cojo_annotate.R), used to annotate SNP rows
##     with snp_easy
## Outputs:
##   data_prep/study_cohort_genetics/output/hla_alleles_snps_regression.csv - association
##     summary statistics, read by figure4_suppfigure1.R and table2_supptable1.R
##     (data_prep/nbb_replication_cohort/3_microglia_neuropathology_regression.R does not use this file - it fits its own
##     microglia/neuropathology regressions directly on the individual-level data)
## =============================================================================

library(data.table)

# Run from the repository root, e.g. `Rscript data_prep/study_cohort_genetics/7_allele_snp_regression.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# glm() can silently drop the `dummy` (allele) term for a group with too
# little variation in it (e.g. a rare allele), in which case it's absent
# from the coefficient table; na_if_empty() converts that empty extraction
# to NA rather than letting it corrupt row lengths.
na_if_empty <- function(x) if (length(x) == 0) NA_real_ else x

# Fits one logistic regression per allele/SNP dosage column in `dosages`
# (dosage ~ dummy + PC1-5), under an additive genetic model. Returns one row
# per allele/SNP with beta, SE, OR (95% CI), p-value, case/control counts and
# allele frequencies, labelled with `group_name`.
run_allele_regression <- function(dosages, group_name) {
    coefficients_df <- data.frame()
    allele_cols <- colnames(dosages)[grepl("^chr6:|\\*|^DR[0-9]$", colnames(dosages))]
    for (allele in allele_cols) {
        dosages_allele <- dosages[!is.na(dosages[[allele]]), ]
        dosages_allele$dummy <- dosages_allele[[allele]]
        model <- try(glm(phenotype ~ dummy + PC1 + PC2 + PC3 + PC4 + PC5, data = dosages_allele, family = 'binomial'), silent = TRUE)
        if (inherits(model, "try-error")) {
            next
        }
        coefficients <- coef(summary(model))
        dummy_coefficient <- na_if_empty(coefficients[grep("dummy", rownames(coefficients)), "Estimate"])
        dummy_standard_error <- na_if_empty(coefficients[grep("dummy", rownames(coefficients)), "Std. Error"])
        dummy_p_value <- na_if_empty(coefficients[grep("dummy", rownames(coefficients)), "Pr(>|z|)"])
        odds_ratio <- exp(dummy_coefficient)
        lower_ci <- exp(dummy_coefficient - 1.96 * dummy_standard_error)
        upper_ci <- exp(dummy_coefficient + 1.96 * dummy_standard_error)
        num_cases <- sum(dosages_allele$phenotype == 1)
        num_controls <- sum(dosages_allele$phenotype == 0)
        c_cases <- sum(dosages_allele$phenotype == 1 & dosages_allele$dummy > 0, na.rm = TRUE)
        c_controls <- sum(dosages_allele$phenotype == 0 & dosages_allele$dummy > 0, na.rm = TRUE)
        homo_cases <- sum(dosages_allele$phenotype == 1 & dosages_allele$dummy == 2, na.rm = TRUE)
        homo_controls <- sum(dosages_allele$phenotype == 0 & dosages_allele$dummy == 2, na.rm = TRUE)
        het_cases <- sum(dosages_allele$phenotype == 1 & dosages_allele$dummy == 1, na.rm = TRUE)
        het_controls <- sum(dosages_allele$phenotype == 0 & dosages_allele$dummy == 1, na.rm = TRUE)
        count_cases <- 2 * homo_cases + het_cases
        count_controls <- 2 * homo_controls + het_controls
        allele_info <- data.frame(
            allele = allele,
            beta = dummy_coefficient,
            standard_error = dummy_standard_error,
            odds_ratio = odds_ratio,
            lower_ci = lower_ci,
            upper_ci = upper_ci,
            p_value = dummy_p_value,
            num_cases = num_cases,
            num_controls = num_controls,
            c_cases = c_cases,
            c_controls = c_controls,
            homo_cases = homo_cases,
            homo_controls = homo_controls,
            het_cases = het_cases,
            het_controls = het_controls,
            count_cases = count_cases,
            count_controls = count_controls,
            freq_cases = count_cases / (2 * num_cases),
            freq_controls = count_controls / (2 * num_controls),
            freq_total = (count_cases + count_controls) / (2 * (num_cases + num_controls)),
            association = group_name
        )
        coefficients_df <- rbind(coefficients_df, allele_info)
    }
    coefficients_df
}

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

CONTROL_DIAGNOSES <- c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD", "Control_MS")
AD_DIAGNOSES <- c("Probable_AD", "AD_path", "Possible AD")

# Per-participant HLA-II allele/SNP dosages (one-field and two-field
# resolution), with diagnosis, sex and age
alleles_snps <- fread('data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv', h = T, sep = ';', stringsAsFactors = F)

diagnosis_counts <- c(
    Centenarian = sum(alleles_snps$Diagnosis == "Centenarian", na.rm = TRUE),
    Controls = sum(alleles_snps$Diagnosis %in% CONTROL_DIAGNOSES, na.rm = TRUE),
    AD = sum(alleles_snps$Diagnosis %in% AD_DIAGNOSES, na.rm = TRUE)
)
message("Cohort composition - Centenarians: ", diagnosis_counts["Centenarian"],
        ", Controls: ", diagnosis_counts["Controls"],
        ", AD: ", diagnosis_counts["AD"],
        ", Total: ", sum(diagnosis_counts))

# Median age at blood collection (±SD) and % female per group, as reported
# in the manuscript.
report_group_demographics <- function(data, group_name, diagnosis_filter) {
    g <- data[data$Diagnosis %in% diagnosis_filter, ]
    age <- as.numeric(g$age)
    pct_female <- 100 * mean(g$sex == "F", na.rm = TRUE)
    message(sprintf(
        "%s: median age at blood collection %.1f±%.1f (N=%d, %.1f%% female)",
        group_name, median(age, na.rm = TRUE), sd(age, na.rm = TRUE),
        nrow(g), pct_female
    ))
}
report_group_demographics(alleles_snps, "CHCs", "Centenarian")
report_group_demographics(alleles_snps, "Controls", CONTROL_DIAGNOSES)
report_group_demographics(alleles_snps, "AD patients", AD_DIAGNOSES)

## -----------------------------------------------------------------------------
## 2. Controls vs. AD
## -----------------------------------------------------------------------------

reg_ctrl_ad <- subset(alleles_snps, Diagnosis %in% c(CONTROL_DIAGNOSES, AD_DIAGNOSES))
reg_ctrl_ad$phenotype <- ifelse(reg_ctrl_ad$Diagnosis %in% CONTROL_DIAGNOSES, 0, 1)
coeff_ctrl_ad <- run_allele_regression(reg_ctrl_ad, "ctrl_ad")

## -----------------------------------------------------------------------------
## 3. Centenarians (CHC) vs. controls
## -----------------------------------------------------------------------------

reg_cent_ctrl <- subset(alleles_snps, Diagnosis %in% c("Centenarian", CONTROL_DIAGNOSES))
reg_cent_ctrl$phenotype <- ifelse(reg_cent_ctrl$Diagnosis == "Centenarian", 0, 1)
coeff_cent_ctrl <- run_allele_regression(reg_cent_ctrl, "cent_ctrl")

## -----------------------------------------------------------------------------
## 4. Centenarians (CHC) vs. AD
## -----------------------------------------------------------------------------

reg_ad_cent <- subset(alleles_snps, Diagnosis %in% c("Centenarian", AD_DIAGNOSES))
reg_ad_cent$phenotype <- ifelse(reg_ad_cent$Diagnosis == "Centenarian", 0, 1)
coeff_ad_cent <- run_allele_regression(reg_ad_cent, "cent_ad")

## -----------------------------------------------------------------------------
## 5. Combine all results and annotate with COJO haplotype labels
## -----------------------------------------------------------------------------

all_coefficients <- rbind(coeff_ctrl_ad, coeff_cent_ctrl, coeff_ad_cent)

# Add the easy-to-read snp_easy label (e.g. "rs35472547") for COJO haplotype
# lead SNPs; NA for all other alleles
cojo_haps_genes <- read.csv2("data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv")
all_coefficients <- merge(all_coefficients, cojo_haps_genes[, c("snp_dos", "snp_easy")], by.x = "allele", by.y = "snp_dos", all.x = TRUE)

## -----------------------------------------------------------------------------
## 6. Male-only DRB1*15:01 frequency and association, for the Discussion's
##    replication check against Mendoza-Mejía et al.'s male-specific finding
## -----------------------------------------------------------------------------

# Reference calculation only (not part of the main association analysis, which
# is not sex-stratified) - the combined-sex frequency is verified against
# Supplementary Table 1 (Hap-Y DRB1*15:01 Freq.CHC/Freq.CNTR).
freq <- function(x) sum(x, na.rm = TRUE) / (2 * sum(!is.na(x)))
chc_male <- alleles_snps$Diagnosis == "Centenarian" & alleles_snps$sex == "M"
cntr_male <- alleles_snps$Diagnosis %in% CONTROL_DIAGNOSES & alleles_snps$sex == "M"
message(sprintf(
    "DRB1*15:01 frequency, males only: CHC=%.4f (N=%d), Controls=%.4f (N=%d)",
    freq(alleles_snps[["DRB1*15:01"]][chc_male]), sum(chc_male),
    freq(alleles_snps[["DRB1*15:01"]][cntr_male]), sum(cntr_male)
))

# Association test for the same comparison, using the main analysis's model
# (dosage ~ group + PC1-5, additive logistic)
sub_male <- alleles_snps[(alleles_snps$Diagnosis == "Centenarian" | alleles_snps$Diagnosis %in% CONTROL_DIAGNOSES) & alleles_snps$sex == "M", ]
sub_male$phenotype <- ifelse(sub_male$Diagnosis == "Centenarian", 0, 1)
sub_male$dummy <- sub_male[["DRB1*15:01"]]
model_male <- glm(phenotype ~ dummy + PC1 + PC2 + PC3 + PC4 + PC5, data = sub_male, family = "binomial")
co_male <- coef(summary(model_male))
message(sprintf("DRB1*15:01, male CHC vs male Control: beta=%.4f, SE=%.4f, P=%.4g, OR=%.3f",
                 co_male["dummy", "Estimate"], co_male["dummy", "Std. Error"],
                 co_male["dummy", "Pr(>|z|)"], exp(co_male["dummy", "Estimate"])))

write.csv2(all_coefficients, "data_prep/study_cohort_genetics/output/hla_alleles_snps_regression.csv", row.names = FALSE)
