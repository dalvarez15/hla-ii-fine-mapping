## =============================================================================
## Haplotype LD (Figure 3): pairwise LD matrix between HLA-II alleles and
## AD/longevity SNPs
## =============================================================================
##
## Reads private data not included in this repository - see DATA_ACCESS.md.
##
## Inputs:
##   data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv - combined
##     SNP + allele dosages (from 5_get_snps_alleles_dosages.R). That file also carries one-field
##     allele columns (e.g. "DRB1*04") for 7_allele_snp_regression.R; the column selection below
##     picks out only SNPs, two-field alleles and DR-broad codes for the LD matrix, excluding those.
## Outputs:
##   data_prep/study_cohort_genetics/output/haplotype_ld/hla_alleles_snps_ld_matrix.{ld,bim} -
##     LD matrix and matching allele/SNP IDs, read directly by figure3.R
## =============================================================================

library(data.table)

# Run from the repository root, e.g. `Rscript data_prep/study_cohort_genetics/6_haplotype_ld.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

OUTDIR <- "data_prep/study_cohort_genetics/output/haplotype_ld"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. Load combined dosages and report cohort composition
## -----------------------------------------------------------------------------

alleles_snps <- fread("data_prep/study_cohort_genetics/output/hla_alleles_snps_dosages_study_cohort.csv", h = T, sep = ';', stringsAsFactors = F)

message("Centenarians: ", sum(alleles_snps$Diagnosis == "Centenarian"),
        ", Controls: ", sum(alleles_snps$Diagnosis %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD", "Control_MS")),
        ", AD cases: ", sum(alleles_snps$Diagnosis %in% c("Probable_AD", "AD_path", "Possible AD")))

## -----------------------------------------------------------------------------
## 2. Build a PLINK dosage (.traw) file: one row per allele/SNP
## -----------------------------------------------------------------------------

# All SNP (chr6:...), two-field HLA allele (Locus*NN:NN) and DR-broad
# (DR1/DR2/...) dosage columns - excludes the one-field allele columns (e.g.
# "DRB1*04") also present in this file, since Figure 3 is two-field only
dosage_matrix <- alleles_snps[, c("ID_GWAS",
                                   grep("^chr6:", colnames(alleles_snps), value = TRUE),
                                   grep("\\*[0-9]+:[0-9]+$|^DR[0-9]$", colnames(alleles_snps), value = TRUE)),
                               with = FALSE]
sample_names <- dosage_matrix$ID_GWAS
dosage_matrix$ID_GWAS <- NULL

dosage_matrix_t <- data.frame(t(dosage_matrix))
colnames(dosage_matrix_t) <- paste("0", sample_names, sep = "_")
dosage_matrix_t$SNP <- rownames(dosage_matrix_t)
dosage_matrix_t$A1 <- "P"
dosage_matrix_t$A2 <- "A"
# Move SNP/A1/A2 to the front, as required by PLINK's .traw format
last_three_indices <- (ncol(dosage_matrix_t) - 2):ncol(dosage_matrix_t)
remaining_indices <- 1:(ncol(dosage_matrix_t) - 3)
df_reordered <- dosage_matrix_t[, c(last_three_indices, remaining_indices)]

BASE <- file.path(OUTDIR, "hla_alleles_snps_ld_matrix")
write.table(df_reordered, paste0(BASE, ".traw"), quote = F, row.names = F, sep = "\t")

fam <- data.frame(FID = 0, IID = sample_names, PAT = 0, MAT = 0, SEX = 0, PHENO = 0)
write.table(fam, paste0(BASE, ".tfam"), quote = F, row.names = F, col.names = FALSE, sep = "\t")

## -----------------------------------------------------------------------------
## 3. Convert to PLINK binary format and compute the LD matrix
## -----------------------------------------------------------------------------

system(sprintf('%s --import-dosage %s.traw id-delim=_ format=1 --fam %s.tfam --real-ref-alleles --make-pgen --out %s', PLINK2_PATH, BASE, BASE, BASE))
system(sprintf('%s --pfile %s --make-bed --out %s', PLINK2_PATH, BASE, BASE))
system(sprintf('%s --bfile %s --r2 square --out %s', PLINK1_PATH, BASE, BASE))

# Keep only the final .ld + .bim; drop the PLINK working files
file.remove(Sys.glob(paste0(BASE, ".", c("traw", "tfam", "pgen", "psam", "pvar", "nosex", "bed", "fam", "log"))))
