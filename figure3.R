## =============================================================================
## Figure 3: LD between the COJO haplotype lead SNPs and previously reported
## AD/longevity GWAS SNPs (3A), and between the three haplotypes that overlap
## the literature (Hap-B, Hap-R, Hap-Y) and the imputed two-field HLA class II
## alleles they tag (3B)
## =============================================================================
##
## Inputs:
##   data_prep/study_cohort_genetics/output/haplotype_ld/hla_alleles_snps_ld_matrix.{ld,bim} -
##     PLINK LD (r2) matrix and identifiers for HLA-region SNPs and imputed HLA alleles;
##     from 6_haplotype_ld.R
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - COJO haplotype
##     lead SNPs (minor-allele oriented) with easy-to-read SNP labels; from 3_cojo_annotate.R
## Outputs:
##   figures/figure3A_all_haplotypes.pdf
##   figures/figure3A_R.pdf
##   figures/figure3B_all_alleles.pdf
##   figures/figure3B_R.pdf
## =============================================================================

library(data.table)
library(corrplot)

# Run from the repository root, e.g. `Rscript figure3.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

# PLINK LD (r2) matrix for HLA-region SNPs and imputed HLA alleles
ld_matrix <- fread('data_prep/study_cohort_genetics/output/haplotype_ld/hla_alleles_snps_ld_matrix.ld', h = F)
ld_bim <- fread('data_prep/study_cohort_genetics/output/haplotype_ld/hla_alleles_snps_ld_matrix.bim', h = F)
ld_ids <- ld_bim$V2 # second BIM column holds the SNP/allele identifiers
ld_matrix <- as.matrix(ld_matrix)
rownames(ld_matrix) <- ld_ids
colnames(ld_matrix) <- ld_ids

# COJO haplotype lead SNPs, with easy-to-read SNP labels (snp_easy)
cojo_haplotypes <- read.csv2("data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv")

## -----------------------------------------------------------------------------
## 2. Annotate SNP/allele identifiers with haplotype and literature labels
## -----------------------------------------------------------------------------

# Reshape the LD matrix to long format, keeping only cross-gene pairs (the
# HLA-II genes are each imputed at high resolution, so within-gene LD between
# alleles of the same gene is not informative here)
ld_long <- as.data.frame(as.table(ld_matrix))
colnames(ld_long) <- c("allele1", "allele2", "r2")
# DR-broad codes (DR1/DR2/...) have no "*" to split on but still belong to
# the HLA-DRB1 locus, same as its two-field alleles (DRB1*01:01 etc.)
normalize_locus <- function(x) ifelse(grepl("^DR[0-9]$", x), "DRB1", sub("\\*.*", "", x))
ld_long$gene1 <- normalize_locus(ld_long$allele1)
ld_long$gene2 <- normalize_locus(ld_long$allele2)
ld_long <- ld_long[ld_long$gene1 != ld_long$gene2, ]
ld_long <- ld_long[order(-ld_long$r2), ]

# Append the COJO haplotype's easy-to-read SNP label to each SNP identifier
ld_long$allele1_ann <- sub("_[^_]+$", "", ld_long$allele1)
ld_long$allele2_ann <- sub("_[^_]+$", "", ld_long$allele2)
ld_long <- merge(ld_long, cojo_haplotypes[, c("SNP", "snp_easy")], by.x = "allele1_ann", by.y = "SNP", all.x = TRUE)
colnames(ld_long)[colnames(ld_long) == "snp_easy"] <- "ann1"
ld_long <- merge(ld_long, cojo_haplotypes[, c("SNP", "snp_easy")], by.x = "allele2_ann", by.y = "SNP", all.x = TRUE)
colnames(ld_long)[colnames(ld_long) == "snp_easy"] <- "ann2"

ld_long$allele1_ann <- ifelse(is.na(ld_long$ann1), ld_long$allele1_ann, paste0(ld_long$allele1_ann, " / ", ld_long$ann1))
ld_long$allele2_ann <- ifelse(is.na(ld_long$ann2), ld_long$allele2_ann, paste0(ld_long$allele2_ann, " / ", ld_long$ann2))

# Append literature annotations for SNPs previously reported in AD/longevity
# GWAS, so the correlation matrix rows/columns display both the SNP ID and
# its source study
annotate_previous_study <- function(x, snp_id, label) {
  ifelse(x == snp_id, paste0(x, " / ", label), x)
}
for (col in c("allele1_ann", "allele2_ann")) {
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32615322:A:G", "Bellenguez AD GWAS 2022")
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32610753:C:A", "Lambert AD GWAS 2013")
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32605638:A:G", "LeGuen ADPD GWAS 2023")
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32622991:C:A", "Joshi Longevity GWAS 2017")
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32684419:T:C", "Timmers Parental-Lifespan GWAS 2019")
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32607629:A:T", "Kunkle AD GWAS 2019")
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32411770:C:T", "Thomassen APOE GWAS 2026")
  ld_long[[col]] <- annotate_previous_study(ld_long[[col]], "chr6:32464090:G:T", "Thomassen APOE GWAS 2026")
}

# Collapse back to a square matrix of annotated allele/SNP labels
cor_matrix <- with(ld_long, tapply(r2, list(allele1_ann, allele2_ann), function(values) values[1]))
cor_matrix[is.na(cor_matrix)] <- 0 # NAs as 0s to avoid corrplot errors

## -----------------------------------------------------------------------------
## 3. Figure 3A: LD between COJO haplotype lead SNPs and literature GWAS SNPs
## -----------------------------------------------------------------------------

ld_colors <- c("#313695", "#74ADD1", "#4DAC26", "#FDAE61", "#D73027")  # Dark blue, light blue, light green, orange, red
color_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)  # Breakpoints for colors
col_fun <- colorRampPalette(ld_colors)

snp_rows <- grepl("^chr6", rownames(cor_matrix))
snp_cols <- grepl("^chr6", colnames(cor_matrix))
snp_rows_cojo <- snp_rows & !grepl("GWAS", rownames(cor_matrix))
snp_cols_gwas <- snp_cols & grepl("GWAS", colnames(cor_matrix))
cor_matrix_snps_gwas <- cor_matrix[snp_rows_cojo, snp_cols_gwas]

desired_cols_order <- c(
  "chr6:32622991:C:A / Joshi Longevity GWAS 2017",
  "chr6:32684419:T:C / Timmers Parental-Lifespan GWAS 2019",
  "chr6:32615322:A:G / Bellenguez AD GWAS 2022",
  "chr6:32605638:A:G / LeGuen ADPD GWAS 2023",
  "chr6:32411770:C:T / Thomassen APOE GWAS 2026",
  "chr6:32610753:C:A / Lambert AD GWAS 2013",
  "chr6:32607629:A:T / Kunkle AD GWAS 2019",
  "chr6:32464090:G:T / Thomassen APOE GWAS 2026"
)
desired_rows_order <- c(
  "chr6:32064966:G:A / Hap-1",
  "chr6:32592593:G:T / Hap-B",
  "chr6:32447376:C:T / Hap-R",
  "chr6:32609869:A:G / Hap-Y",
  "chr6:32707290:T:C / Hap-5",
  "chr6:32772190:A:T / Hap-6"
)

# All six COJO haplotype lead SNPs against the literature GWAS SNPs
pdf("figures/figure3A_all_haplotypes.pdf")
corrplot(as.matrix(cor_matrix_snps_gwas[desired_rows_order, desired_cols_order]), method = "color", is.corr = FALSE, tl.cex = 0.6,
         col = col_fun(length(color_breaks) - 1), tl.col = "black",
         addgrid.col = "white", col.lim = c(0, 1), na.label = " ", addCoef.col = "black", number.cex = 0.5)
dev.off()

# Only the three haplotype lead SNPs with at least one correlation > 0.2 with
# a literature GWAS SNP (Hap-B, Hap-R, Hap-Y)
haplotypes_in_literature <- c(
  "chr6:32592593:G:T / Hap-B",
  "chr6:32447376:C:T / Hap-R",
  "chr6:32609869:A:G / Hap-Y"
)
cor_matrix_snps_gwas_pres <- cor_matrix_snps_gwas[haplotypes_in_literature, desired_cols_order, drop = FALSE]


pdf("figures/figure3A_R.pdf")
corrplot(as.matrix(cor_matrix_snps_gwas_pres), method = "color", is.corr = FALSE, tl.cex = 0.6,
         col = col_fun(length(color_breaks) - 1), tl.col = "black",
         addgrid.col = "white", col.lim = c(0, 1), na.label = " ", addCoef.col = "black", number.cex = 0.5)
dev.off()

## -----------------------------------------------------------------------------
## 4. Figure 3B: LD between COJO haplotype lead SNPs and imputed HLA alleles
## -----------------------------------------------------------------------------

# Restrict to two-field HLA class II alleles
class_ii_two_field <- grepl(
  "^(DRA|DRB[1-5]|DQA[1-2]|DQB[1-2]|DPA[1-2]|DPB[1-2])\\*[0-9]+:[0-9]+$",
  colnames(cor_matrix)
)
cor_matrix_cojo_snps_alleles <- cor_matrix[haplotypes_in_literature, class_ii_two_field, drop = FALSE]

# Keep only alleles with at least one correlation > 0.2 with a haplotype lead SNP
cols_to_keep <- apply(cor_matrix_cojo_snps_alleles, 2, function(col) any(col > 0.2))
cor_matrix_cojo_snps_alleles <- cor_matrix_cojo_snps_alleles[, cols_to_keep]

pdf("figures/figure3B_all_alleles.pdf")
corrplot(as.matrix(cor_matrix_cojo_snps_alleles), method = "color", is.corr = FALSE, tl.cex = 0.6,
         col = col_fun(length(color_breaks) - 1), tl.col = "black",
         addgrid.col = "white", col.lim = c(0, 1), na.label = " ", addCoef.col = "black", number.cex = 0.5)
dev.off()

# Restrict to one allele per gene per haplotype (excluding second-field
# alleles that are in LD with an already-selected allele of the same
# haplotype, i.e., DQA1*03:03 with Hap-B's DQA1*03:01, r2=0.33 < 0.60; and
# DRB1*04:04 with Hap-B's DRB1*04:01, r2=0.20 < 0.62), grouped by haplotype
desired_col_order <- c(
  "DRB1*04:01", "DQA1*03:01", "DQB1*03:02", # Hap-B
  "DRB1*01:01", "DQA1*01:01", "DQB1*05:01", # Hap-R
  "DRB1*15:01", "DQA1*01:02", "DQB1*06:02"  # Hap-Y
)
cor_matrix_cojo_snps_alleles <- cor_matrix_cojo_snps_alleles[, desired_col_order, drop = FALSE]

pdf("figures/figure3B_R.pdf")
corrplot(as.matrix(cor_matrix_cojo_snps_alleles), method = "color", is.corr = FALSE, tl.cex = 0.6,
         col = col_fun(length(color_breaks) - 1), tl.col = "black",
         addgrid.col = "white", col.lim = c(0, 1), na.label = " ", addCoef.col = "black", number.cex = 0.5)
dev.off()
