## =============================================================================
## Annotate the r2=0.8 COJO haplotypes with gene context, external AD-GWAS
## overlap, the manuscript's Hap-1/Hap-R/Hap-B/Hap-Y/Hap-5/Hap-6 label per
## haplotype (Table 1), and re-orient every SNP to its minor allele for
## consistent LD/dosage column naming downstream.
## =============================================================================
##
## Inputs:
##   data_prep/study_cohort_genetics/output/cojo/summary_cojo.txt - COJO results (from 1_clumping_cojo.R)
##   data_prep/reference_data/ensembl_gene_coordinates_chr6_mhc.txt - for gene-context annotation
##   AD_MEGA_META_HITS (raw_input_data/data_paths.R) - external AD-GWAS meta-analysis
##     index-variant hits, for GWAS-overlap annotation (aggregate, not individual-level)
## Outputs:
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - annotated COJO
##     haplotypes, re-oriented to the minor allele; read directly by figure1_figure2.R, figure3.R,
##     table1.R, 7_allele_snp_regression.R and ../nbb_replication_cohort/3_microglia_neuropathology_regression.R
## =============================================================================

library(data.table)
library(dplyr)

# Run from the repository root, e.g. `Rscript data_prep/study_cohort_genetics/3_cojo_annotate.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

source("raw_input_data/data_paths.R")

## -----------------------------------------------------------------------------
## 1. Load the COJO results
## -----------------------------------------------------------------------------

cojo_haplotypes <- fread("data_prep/study_cohort_genetics/output/cojo/summary_cojo.txt", data.table = FALSE)
cojo_haplotypes$SNP <- paste0("chr", cojo_haplotypes$SNP)

## -----------------------------------------------------------------------------
## 2. Annotate with gene context and external AD-GWAS overlap
## -----------------------------------------------------------------------------

ensembl_genes <- read.delim("data_prep/reference_data/ensembl_gene_coordinates_chr6_mhc.txt", header = TRUE, sep = "\t")
colnames(ensembl_genes) <- c("gene_id", "chromosome", "start", "end", "gene_name")
ensembl_genes <- ensembl_genes %>%
  filter(!is.na(gene_name) & gene_name != "")

# For each lead SNP, report the overlapping gene(s), or the nearest genes on
# either side if it falls between genes
cojo_haplotypes <- cojo_haplotypes %>%
    rowwise() %>%
    mutate(gene = {
        gene_match <- ensembl_genes %>%
            filter(chromosome == Chr, bp >= start, bp <= end) %>%
            pull(gene_name)
        if (length(gene_match) > 1) {
            paste(gene_match, collapse = ", ")
        } else if (length(gene_match) == 1) {
            gene_match[1]
        } else {
            before_gene <- ensembl_genes %>%
                filter(chromosome == Chr, start < bp) %>%
                arrange(desc(start)) %>%
                slice(1) %>%
                pull(gene_name)
            after_gene <- ensembl_genes %>%
                filter(chromosome == Chr, start > bp) %>%
                arrange(start) %>%
                slice(1) %>%
                pull(gene_name)
            paste(before_gene, after_gene, sep = ", ")
        }
    }) %>%
    ungroup()
cojo_haplotypes$cojo <- TRUE

# Flag lead SNPs that are also reported index variants in an external AD-GWAS
# meta-analysis, and carry over that analysis's joint effect estimate
gwas_hits <- fread(cmd = paste("zcat", AD_MEGA_META_HITS), header = TRUE, sep = "\t", data.table = FALSE)
cojo_haplotypes$gwas <- cojo_haplotypes$SNP %in% unique(gwas_hits$index_variant)
cojo_haplotypes <- merge(cojo_haplotypes, gwas_hits[, c("index_variant", "b_cojo_slct")], by.x = "SNP", by.y = "index_variant", all.x = TRUE)

## -----------------------------------------------------------------------------
## 3. Assign the manuscript's Hap-X label to each of the six candidate haplotypes
## -----------------------------------------------------------------------------

cojo_haplotypes$snp_easy <- NA
cojo_haplotypes[cojo_haplotypes$SNP == "chr6:32064966:G:A", "snp_easy"] <- "Hap-1"
cojo_haplotypes[cojo_haplotypes$SNP == "chr6:32447376:C:T", "snp_easy"] <- "Hap-R" # EADB-GWAS-2026 risk lead SNP (rs9469112)
cojo_haplotypes[cojo_haplotypes$SNP == "chr6:32592593:G:T", "snp_easy"] <- "Hap-B" # EADB-GWAS-2026 protective lead SNP (rs35472547)
cojo_haplotypes[cojo_haplotypes$SNP == "chr6:32609869:A:G", "snp_easy"] <- "Hap-Y"
cojo_haplotypes[cojo_haplotypes$SNP == "chr6:32707290:T:C", "snp_easy"] <- "Hap-5"
cojo_haplotypes[cojo_haplotypes$SNP == "chr6:32772190:A:T", "snp_easy"] <- "Hap-6"

## -----------------------------------------------------------------------------
## 4. Build the dosage-column identifier (reference-allele oriented)
## -----------------------------------------------------------------------------

cojo_haplotypes$snp_dos <- paste0(cojo_haplotypes$SNP, "_", cojo_haplotypes$refA)
cojo_haplotypes <- cojo_haplotypes %>%
    mutate(
        altA = sapply(seq_len(n()), function(i) {
            alleles <- strsplit(SNP[i], ":")[[1]]
            allele1 <- alleles[3]
            allele2 <- alleles[4]
            if (refA[i] == allele1) allele2 else if (refA[i] == allele2) allele1 else NA
        })
    )

## -----------------------------------------------------------------------------
## 5. Re-orient to the minor allele and write the final file
## -----------------------------------------------------------------------------

cojo_haplotypes_flipped <- cojo_haplotypes %>%
    mutate(
        flip = freq > 0.5,
        refA_tmp = ifelse(flip, altA, refA),
        altA_tmp = ifelse(flip, refA, altA),
        freq = ifelse(flip, 1 - freq, freq),
        freq_geno = ifelse(flip, 1 - freq_geno, freq_geno),
        b = ifelse(flip, -b, b),
        bJ = ifelse(flip, -bJ, bJ),
        b_cojo_slct = ifelse(flip, -b_cojo_slct, b_cojo_slct),
        refA = refA_tmp,
        altA = altA_tmp
    ) %>%
    select(-refA_tmp, -altA_tmp, -flip)
cojo_haplotypes_flipped$snp_dos <- paste0(cojo_haplotypes_flipped$SNP, "_", cojo_haplotypes_flipped$refA)

write.csv2(cojo_haplotypes_flipped, file = "data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv", row.names = FALSE)
