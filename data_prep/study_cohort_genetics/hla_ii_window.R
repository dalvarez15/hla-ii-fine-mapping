## =============================================================================
## HLA-II analysis window on chromosome 6: the EADB-GWAS-2026 SNPs within
## 0.5 Mb of the locus's two lead SNPs (risk rs9469112, protective
## rs35472547), from 32 Mb onward. With the EADB-GWAS-2026 summary statistics
## this is chr6:32,037,271-33,092,341 (~1 Mb, spanning HLA-DRB1, HLA-DQA1 and
## HLA-DQB1).
## =============================================================================
##
## Sourced by 1_clumping_cojo.R (COJO candidate clumps), 2_regional_ld.R
## (regional LD) and ../../figure1_figure2.R (Figures 1 and 2), so all three
## use the same window.
## =============================================================================

HLA_II_MIN_POS <- 32000000
LEAD_SNP_POS   <- c(32447376, 32592593) # rs9469112 (risk), rs35472547 (protective)
LEAD_FLANK     <- 500000

# Returns the window's first and last SNP position, given the positions of
# the chr6 SNPs in the summary statistics
hla_ii_window <- function(pos) {
    pos <- pos[pos >= min(LEAD_SNP_POS) - LEAD_FLANK &
               pos <= max(LEAD_SNP_POS) + LEAD_FLANK &
               pos >= HLA_II_MIN_POS]
    c(start = min(pos), end = max(pos))
}
