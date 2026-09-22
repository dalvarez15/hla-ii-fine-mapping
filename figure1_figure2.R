## =============================================================================
## Figure 1: regional association plot for the HLA/MHC class II region,
## colored/shaped by LD (r2) with the two EADB-GWAS-2026 top SNPs (protective:
## rs35472547; risk: rs9469112), with previously reported AD/longevity SNPs
## highlighted.
## Figure 2: A) all suggestive (P<P_SUGGESTIVE) clumping haplotypes; B) the
## six COJO candidate haplotypes and their lead SNPs; C) each haplotype
## lead SNP's FDR-adjusted joint (COJO) p-value. Also reports the
## region/clump/SNP counts and percentages cited in the Results text and
## Figure 2 legend.
## =============================================================================
##
## Inputs:
##   data_prep/reference_data/gwas_summary_stats_eadb2026.txt - EADB-GWAS-2026 summary
##     statistics, aligned with PLINK, trimmed to chr6:32,000,000-34,000,000 (the search window
##     this script's own region-boundary computation below needs) - see DATA_ACCESS.md for why
##     only this narrow window is committed here rather than the full chr6 (or genome-wide) file
##   data_prep/study_cohort_genetics/output/clumping_r2_080.clumps - from 1_clumping_cojo.R
##   data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv - from 3_cojo_annotate.R
##   data_prep/reference_data/ensembl_gene_coordinates_chr6_mhc.txt - for the HLA class II gene track
##   data_prep/study_cohort_genetics/output/regional_ld/regional_ld_pairs.ld - from 2_regional_ld.R
## Outputs:
##   figures/figure1_R.pdf
##   figures/figure2_R.pdf
## =============================================================================

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Run from the repository root, e.g. `Rscript figure1_figure2.R` from the top-level directory,
# or uncomment and edit the line below to set the working directory explicitly.
# setwd("/path/to/repository")

## -----------------------------------------------------------------------------
## Global parameters
## -----------------------------------------------------------------------------

CHR          <- 6
P_GW         <- 5e-8
# "Suggestive" significance threshold, matching the one used upstream when
# generating the clumps themselves (see data_prep/study_cohort_genetics/1_clumping_cojo.R)
# to flag clumps eligible for COJO follow-up. Used consistently here for the
# SNP-level plot lines/counts and for Figure 2 panel A's clump count/coloring.
P_SUGGESTIVE <- 1e-5

PROTECTIVE_TOP_SNP <- "6:32592593:G:T" # EADB-GWAS-2026 top SNP, rs35472547
RISK_TOP_SNP       <- "6:32447376:C:T" # EADB-GWAS-2026 top SNP, rs9469112

## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------

# Expands each clump in `clumps_df` (one row per clump, with lead SNP `ID`
# and comma-separated secondary members `SP2`) into per-row membership and
# color in `data`, matched by SNP identifier (`data$ID`). Rows not in any
# clump keep `default_clump`/`default_color`.
assign_clump_membership <- function(data, clumps_df, color_for_lead, default_clump, default_color = "grey") {
  clump <- rep(default_clump, nrow(data))
  color <- rep(default_color, nrow(data))
  for (i in seq_len(nrow(clumps_df))) {
    lead <- clumps_df$ID[i]
    members <- unique(c(lead, unlist(strsplit(clumps_df$SP2[i], ","))))
    members <- members[nzchar(members)]
    idx <- which(data$ID %in% members)
    clump[idx] <- lead
    color[idx] <- color_for_lead(lead)
  }
  data.frame(clump = clump, color = color, stringsAsFactors = FALSE)
}

## -----------------------------------------------------------------------------
## 1. Load input data
## -----------------------------------------------------------------------------

gwas_summary_stats <- fread(
  "data_prep/reference_data/gwas_summary_stats_eadb2026.txt"
)
colnames(gwas_summary_stats)[colnames(gwas_summary_stats) == "#CHROM"] <- "CHROM"

clumps <- fread(
  "data_prep/study_cohort_genetics/output/clumping_r2_080.clumps"
)

cojo_haps <- read.csv2(
  "data_prep/study_cohort_genetics/output/cojo_haplotypes_minor_allele.csv"
)
cojo_haps$ID <- sub("^chr", "", cojo_haps$SNP)
hap_rs <- c("Hap-B" = "rs35472547", "Hap-R" = "rs9469112", "Hap-Y" = "rs4335021")
cojo_haps$rs <- unname(hap_rs[cojo_haps$snp_easy])

ensembl_genes <- fread(
  "data_prep/reference_data/ensembl_gene_coordinates_chr6_mhc.txt",
  col.names = c("gene_id", "chromosome", "start", "end", "gene")
)

# Previously reported SNPs (AD or longevity literature) that fall within the
# region plotted here, used to flag them in Figure 1
previous_studies_lookup <- data.frame(
    ID = c(
        "6:32615322:A:G", "6:32610753:C:A", "6:32605638:A:G", "6:32607629:A:T",
        "6:32622991:C:A", "6:32684419:T:C",
        "6:32411770:C:T", "6:32464090:G:T"
    ),
    rs = c(
        "rs6605556", "rs9271192", "rs601945", "rs9271058",
        "rs34831921", "rs9275152",
        "rs17208902", "rs9268888"
    ),
    study = c(
        "Bellenguez 2022", "Lambert 2013", "Le Guen 2024", "Kunkle 2019",
        "Joshi 2017", "Timmers 2019",
        "Thomassen-1 2026", "Thomassen-2 2026"
    ),
    stringsAsFactors = FALSE
)

## -----------------------------------------------------------------------------
## 2. Preprocess GWAS summary statistics for the HLA/MHC region
## -----------------------------------------------------------------------------

# Region spans the HLA class II genes, bounded by the first SNP at or after
# 32 Mb and the last SNP at or before 34 Mb on chromosome 6
plot_start <- min(gwas_summary_stats$POS[gwas_summary_stats$CHROM == CHR & gwas_summary_stats$POS >= 32000000])
plot_end <- max(gwas_summary_stats$POS[gwas_summary_stats$CHROM == CHR & gwas_summary_stats$POS >= 33000000 & gwas_summary_stats$POS <= 34000000])

fullstats_filtered <- gwas_summary_stats %>%
  filter(CHROM == CHR, POS >= plot_start, POS <= plot_end) %>%
  filter(!duplicated(POS)) # keep one SNP per position

## -----------------------------------------------------------------------------
## 3. Identify COJO-significant clumps and assign clump colors
## -----------------------------------------------------------------------------

clumps_sig <- clumps %>% filter(ID %in% cojo_haps$ID)

# Fixed colors for the three top GWAS SNPs, gradient palette for the rest
special_colors <- c(
  "6:32447376:C:T" = "red",   # RISK_TOP_SNP
  "6:32592593:G:T" = "blue",  # PROTECTIVE_TOP_SNP
  "6:32609869:A:G" = "gold"
)
remaining_ids <- setdiff(clumps_sig$ID, names(special_colors))
palette <- brewer.pal(max(7, length(remaining_ids)), "Dark2")
clump_color_map <- c(
  special_colors,
  setNames(palette[seq_along(remaining_ids)], remaining_ids)
)

clump_membership <- assign_clump_membership(
  fullstats_filtered, clumps_sig,
  color_for_lead = function(lead) clump_color_map[[lead]],
  default_clump = "none"
)
fullstats_filtered$clump <- clump_membership$clump
fullstats_filtered$color <- clump_membership$color

## -----------------------------------------------------------------------------
## 4. Gene track (HLA class II genes, shared by Figures 1 and 2)
## -----------------------------------------------------------------------------

hla_genes <- ensembl_genes %>%
  filter(gene %in% c("HLA-DQA1", "HLA-DQB1", "HLA-DRB1")) %>%
  mutate(y = 0.8, mid = (start + end) / 2, length = end - start)

gene_track <- ggplot(hla_genes) +
    geom_segment(
        aes(x = start, xend = end, y = y, yend = y),
        linewidth = 4, color = "black", lineend = "round"
    ) +
    geom_text(
        aes(x = mid, y = y, label = gene),
        angle = 35,
        hjust = 1,       # pulls the rotated label downward-left
        vjust = 1.2,     # pushes text below the segment
        nudge_y = -0.05, # small downward offset
        size = 4
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0.05, 0.05))
    ) +
    coord_cartesian(xlim = c(plot_start, plot_end)) +
    theme_void() +
    theme(
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )

## -----------------------------------------------------------------------------
## 5. Figure 2 panel B: regional plot colored by COJO-significant clump
## -----------------------------------------------------------------------------

cojo_haps_plot <- ggplot(fullstats_filtered, aes(POS, -log10(P))) +
  geom_hline(yintercept = -log10(P_GW), color = "red") +
  geom_hline(yintercept = -log10(P_SUGGESTIVE), linetype = "dashed") +
  geom_point(aes(color = color), size = 1.6, alpha = 0.85, shape = 16) +
  scale_color_identity() +
    geom_point(
        data = subset(fullstats_filtered, ID %in% cojo_haps$ID & ID != PROTECTIVE_TOP_SNP & ID != RISK_TOP_SNP),
        aes(fill = color),
        shape = 21,
        color = "black",
        size = 3,
        stroke = 0.7,
        show.legend = FALSE
    ) +
    geom_point(
        data = subset(fullstats_filtered, ID %in% PROTECTIVE_TOP_SNP),
        aes(x = POS, y = -log10(P)),
        color = "black",
        fill = "blue",
        shape = 23,
        size = 3,
        stroke = 0.7
    ) +
        geom_point(
        data = subset(fullstats_filtered, ID %in% RISK_TOP_SNP),
        aes(x = POS, y = -log10(P)),
        color = "black",
        fill = "red",
        shape = 24,
        size = 3,
        stroke = 0.7
    ) +
    scale_fill_identity() +
    theme_bw() +
    labs(
        x = "Genomic coordinate (bp), chromosome 6",
        y = "-log10(P GWAS)"
    ) +
    xlim(plot_start, plot_end) +
    theme(axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14))

## -----------------------------------------------------------------------------
## 6. Figure 2 panel C: COJO joint (FDR) p-value of each haplotype lead SNP
## -----------------------------------------------------------------------------

cojo_haps$color <- clump_color_map[cojo_haps$ID]

cojo_haps_pJ_rev <- ggplot(cojo_haps, aes(bp, -log10(pJ_FDR))) +
    geom_hline(yintercept = -log10(0.05), color = "red") +
    geom_point(
        data = subset(cojo_haps, ID != PROTECTIVE_TOP_SNP & ID != RISK_TOP_SNP),
        aes(fill = color),
        color = "black",
        shape = 21,
        size = 3,
        stroke = 0.7,
        show.legend = FALSE
    ) +
    geom_point(
        data = subset(cojo_haps, ID %in% PROTECTIVE_TOP_SNP),
        color = "black",
        fill = "blue",
        shape = 23,
        size = 3,
        stroke = 0.7
    ) +
        geom_point(
        data = subset(cojo_haps, ID %in% RISK_TOP_SNP),
        color = "black",
        fill = "red",
        shape = 24,
        size = 3,
        stroke = 0.7
    ) +
    scale_y_reverse(limits = c(0, 8)) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(plot_start, plot_end)) +
    labs(x = "Genomic coordinate (bp), chromosome 6", y = "-log10(P COJO FDR)") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14))

## -----------------------------------------------------------------------------
## 7. Figure 2 panel A: regional plot colored by all suggestive-significance
##    clumps (P<P_SUGGESTIVE)
## -----------------------------------------------------------------------------

clumps_filtered <- clumps %>% filter(POS >= plot_start, POS <= plot_end, P < P_SUGGESTIVE)
clumps_ids <- unique(clumps_filtered$ID)
stopifnot(all(clumps_ids %in% fullstats_filtered$ID))

# Gradient palette for all suggestive clumps, avoiding colors already used
# for the COJO-significant clumps and skipping near-black/near-white shades
used_colors <- unique(cojo_haps$color)
available_colors <- setdiff(brewer.pal(8, "Set2"), used_colors)
all_palette <- colorRampPalette(available_colors)(length(clumps_ids) + 17)
all_palette <- all_palette[!grepl("#[A-F0-9]{2}[A-F0-9]{2}[0-3][0-9A-F]", all_palette, ignore.case = TRUE)]
if (length(all_palette) < length(clumps_ids)) {
    all_palette <- rep_len(all_palette, length(clumps_ids))
}

clump_membership_all <- assign_clump_membership(
  fullstats_filtered, clumps_filtered,
  color_for_lead = function(lead) all_palette[which(clumps_ids == lead)],
  default_clump = NA_character_
)
fullstats_clumps <- fullstats_filtered
fullstats_clumps$clump_all <- clump_membership_all$clump
fullstats_clumps$color_all <- clump_membership_all$color

# Override with the fixed COJO-significant clump colors, so those clumps
# match panel B
for (i in seq_along(cojo_haps$ID)) {
    idx <- which(fullstats_clumps$clump_all == cojo_haps$ID[i])
    if (length(idx) > 0) {
        fullstats_clumps$color_all[idx] <- cojo_haps$color[i]
    }
}

clumps_plot <- ggplot(fullstats_clumps, aes(POS, -log10(P), color = color_all)) +
    geom_hline(yintercept = -log10(P_GW), color = "red") +
    geom_hline(yintercept = -log10(P_SUGGESTIVE), linetype = "dashed") +
    geom_point(alpha = 0.9, size = 1.5, show.legend = FALSE) +
    scale_color_identity() +
        geom_point(
        data = subset(fullstats_clumps, ID %in% PROTECTIVE_TOP_SNP),
        aes(x = POS, y = -log10(P)),
        color = "black",
        fill = "blue",
        shape = 23,
        size = 3,
        stroke = 0.7
    ) +
        geom_point(
        data = subset(fullstats_clumps, ID %in% RISK_TOP_SNP),
        aes(x = POS, y = -log10(P)),
        color = "black",
        fill = "red",
        shape = 24,
        size = 3,
        stroke = 0.7
    ) +
    scale_fill_identity() +
    labs(
        x = "Genomic coordinate (bp), chromosome 6",
        y = "-log10(P GWAS)"
    ) +
    coord_cartesian(xlim = c(plot_start, plot_end)) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
            axis.text  = element_text(size = 14))

## -----------------------------------------------------------------------------
## 8. Figure 2: combine panels A-C and save
## -----------------------------------------------------------------------------

gene_track_no_x <- gene_track +
    theme(axis.title.x = element_blank(),
                axis.text.x  = element_blank(),
                axis.ticks.x = element_blank())

clumps_plot_no_x <- clumps_plot +
    theme(axis.title.x = element_blank()) +
    annotate("text", x = Inf, y = Inf, label = "A", vjust = 1.5, hjust = 1.5, size = 5, fontface = "bold")

cojo_haps_plot_no_x <- cojo_haps_plot +
    scale_y_reverse() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()) +
    annotate(
        "label",
        x = plot_end - ((plot_end - plot_start) * 0.005),
        y = 9,
        label = paste0(
            "COJO haplotypes lead SNPs:\n\n",
            paste(paste0("   ", cojo_haps$ID), ifelse(is.na(cojo_haps$rs), "", paste0(", ", cojo_haps$rs)), sep = "", collapse = "\n\n")
        ),
        hjust = 1,
        vjust = 1,
        size = 3.6,
        lineheight = 1.2,
        fill = scales::alpha("white", 0.80),
        colour = "black",
        label.padding = grid::unit(0.6, "lines"),
        label.r = grid::unit(0.15, "lines")
    ) +
    annotate("text", x = Inf, y = 0, label = "B", vjust = 1.5, hjust = 1.5, size = 5, fontface = "bold")

cojo_haps_pJ_rev <- cojo_haps_pJ_rev +
    annotate("text", x = Inf, y = 0, label = "C", vjust = 1.5, hjust = 1.5, size = 5, fontface = "bold")

combined_clumps_gene_cojo_pj <- clumps_plot_no_x / gene_track_no_x / cojo_haps_plot_no_x / cojo_haps_pJ_rev +
    plot_layout(heights = c(6, 1, 6, 3))

ggsave("figures/figure2_R.pdf",
       plot = combined_clumps_gene_cojo_pj, width = 12, height = 16, dpi = 300)

## -----------------------------------------------------------------------------
## 9. Figure 1: LD of every SNP to the two EADB-GWAS-2026 top SNPs
## -----------------------------------------------------------------------------

ld_pairs <- fread('data_prep/study_cohort_genetics/output/regional_ld/regional_ld_pairs.ld')
ld_vec_protective <- setNames(ld_pairs$R2[ld_pairs$SNP_A == PROTECTIVE_TOP_SNP], ld_pairs$SNP_B[ld_pairs$SNP_A == PROTECTIVE_TOP_SNP])
ld_vec_risk <- setNames(ld_pairs$R2[ld_pairs$SNP_A == RISK_TOP_SNP], ld_pairs$SNP_B[ld_pairs$SNP_A == RISK_TOP_SNP])
ld_annot_protective <- ld_vec_protective[match(fullstats_filtered$ID, names(ld_vec_protective))]
ld_annot_risk <- ld_vec_risk[match(fullstats_filtered$ID, names(ld_vec_risk))]
fullstats_filtered$ld_r2_top_snp1 <- ld_annot_protective
fullstats_filtered$ld_r2_top_snp2 <- ld_annot_risk
fullstats_filtered$ld_r2_top_snp1_and2 <- pmax(ld_annot_protective, ld_annot_risk)

fullstats_filtered$max_ld_snp <- ifelse(ld_annot_protective > ld_annot_risk, "top_snp1",
                                         ifelse(ld_annot_risk > ld_annot_protective, "top_snp2", "equal"))

# Shape by which top SNP a given SNP is in higher LD with (diamond/triangle),
# or a plain circle when in low LD (<0.2) with both
fullstats_filtered$point_shape <- ifelse(
    fullstats_filtered$ld_r2_top_snp1_and2 < 0.2 | is.na(fullstats_filtered$ld_r2_top_snp1_and2),
    21,  # circle, low LD with both
    ifelse(fullstats_filtered$max_ld_snp == "top_snp1", 23, 24) # diamond: protective, triangle: risk
)

ld_colors <- c("#313695", "#74ADD1", "#4DAC26", "#FDAE61", "#D73027")
color_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

fullstats_filtered$in_previous_studies <- fullstats_filtered$ID %in% previous_studies_lookup$ID

# This SNP's LD was rounded to 0.60 (from 0.599) in Figure 3, so the same
# rounding is applied here for consistency between figures
fullstats_filtered[which(fullstats_filtered$ID == "6:32411770:C:T"), "ld_r2_top_snp1_and2"] <- 0.61 # Thomassen-1

gwas_ld_topsnps <- ggplot(fullstats_filtered, aes(x = POS, y = -log10(P))) +
    geom_hline(yintercept = -log10(P_GW), color = "red", linewidth = 0.8) +
    geom_hline(yintercept = -log10(P_SUGGESTIVE), linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(aes(color = ld_r2_top_snp1_and2, shape = point_shape, fill = ld_r2_top_snp1_and2), size = 1.6, alpha = 0.9) +
    scale_color_stepsn(colors = ld_colors, breaks = color_breaks, limits = c(0, 1), name = expression(LD~(r^2)), guide = guide_colorsteps(show.limits = TRUE)) +
    scale_fill_stepsn(colors = ld_colors, breaks = color_breaks, limits = c(0, 1), name = expression(LD~(r^2))) +
    scale_shape_identity(
        breaks = c(23, 24, 21),
        labels = c(paste0("rs35472547 (", PROTECTIVE_TOP_SNP, ")"),
                 paste0("rs9469112 (", RISK_TOP_SNP, ")"),
                 "Low LD (<0.2) with both top SNPs"),
        name = "LD with EADB-GWAS-2026 top SNP:",
        guide = "legend"
    ) +
    # previous studies keep their shape, but get a black border
    geom_point(
        data = subset(fullstats_filtered, in_previous_studies),
        aes(x = POS, y = -log10(P), shape = point_shape, fill = ld_r2_top_snp1_and2),
        color = "black",
        size = 2.0,
        stroke = 1,
        alpha = 1,
        show.legend = FALSE
    ) +
    geom_point(data = subset(fullstats_filtered, ID %in% PROTECTIVE_TOP_SNP),
               aes(x = POS, y = -log10(P)), color = "black", fill = "blue", shape = 23, size = 3, stroke = 1.5) +
    geom_point(data = subset(fullstats_filtered, ID %in% RISK_TOP_SNP),
               aes(x = POS, y = -log10(P)), color = "black", fill = "red", shape = 24, size = 3, stroke = 1.5) +
    theme_bw() +
    labs(x = "Genomic coordinate (bp), chromosome 6", y = "-log10(P GWAS)") +
    theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.box = "horizontal",
        legend.direction = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
    ) +
    guides(
        color = guide_colorsteps(order = 2, show.limits = TRUE, title.position = "top"),
        fill = "none",
        shape = guide_legend(order = 1, title.position = "top",
                           override.aes = list(
                            fill = "grey70",
                            color = NA,
                            size = 3, alpha = 1))
    )

# y-axis capped at -log10(P)=27.5, with a taller output to match
gwas_ld_topsnps_y30 <- gwas_ld_topsnps +
    coord_cartesian(xlim = c(plot_start, plot_end), ylim = c(0, 27.5)) +
    scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25))
gwas_ld_topsnps_genetrack_y30 <- gwas_ld_topsnps_y30 / gene_track +
  plot_layout(heights = c(7, 1))

ggsave("figures/figure1_R.pdf",
    plot = gwas_ld_topsnps_genetrack_y30, width = 12, height = 8, dpi = 300)

## -----------------------------------------------------------------------------
## 10. Region/clump/SNP counts reported in the Results text and Figure 2 legend
## -----------------------------------------------------------------------------

# Number of clumps in the region reaching each significance threshold, by the
# clump's own (lead SNP) P-value. n_clumps_suggestive reuses clumps_filtered
# (step 7), which is filtered at the same P_SUGGESTIVE threshold.
n_clumps_suggestive <- nrow(clumps_filtered)
n_clumps_gw <- nrow(clumps[clumps$POS >= plot_start & clumps$POS <= plot_end & clumps$P < P_GW, ])
message(sprintf("Clumps within the region at P<%g: %d", P_SUGGESTIVE, n_clumps_suggestive))
message(sprintf("Clumps within the region at P<%g: %d", P_GW, n_clumps_gw))

# Proportion of SNPs in the region in low LD (<0.2) with both top SNPs,
# overall and among suggestive/genome-wide significant SNPs
total_snps <- nrow(fullstats_filtered[!is.na(fullstats_filtered$ld_r2_top_snp1_and2), ])
low_ld_snps <- sum(fullstats_filtered$ld_r2_top_snp1_and2 < 0.2 & !is.na(fullstats_filtered$ld_r2_top_snp1_and2))
suggestive_snps <- sum(fullstats_filtered$P < P_SUGGESTIVE, na.rm = TRUE)
suggestive_low_ld <- sum(fullstats_filtered$P < P_SUGGESTIVE & fullstats_filtered$ld_r2_top_snp1_and2 < 0.2, na.rm = TRUE)
gw_snps <- sum(fullstats_filtered$P < P_GW, na.rm = TRUE)
gw_low_ld <- sum(fullstats_filtered$P < P_GW & fullstats_filtered$ld_r2_top_snp1_and2 < 0.2, na.rm = TRUE)

message(sprintf("Plotted region: chr%d:%d-%d (%.2f Mb); %.0f kb from %s (protective), %.0f kb from %s (risk)",
                CHR, plot_start, plot_end, (plot_end - plot_start) / 1e6,
                (plot_end - as.numeric(str_split_fixed(PROTECTIVE_TOP_SNP, ":", 4)[, 2])) / 1e3, "rs35472547",
                (as.numeric(str_split_fixed(RISK_TOP_SNP, ":", 4)[, 2]) - plot_start) / 1e3, "rs9469112"))
message(sprintf("Total SNPs in region: %d", total_snps))
message(sprintf("SNPs with low LD (r2 < 0.2) to either EADB-GWAS-2026 SNP: %d (%.1f%%)",
                low_ld_snps, 100 * low_ld_snps / total_snps))
message(sprintf("SNPs with suggestive significance (P < %g): %d", P_SUGGESTIVE, suggestive_snps))
message(sprintf("Suggestive SNPs with low LD: %d (%.1f%%)", suggestive_low_ld,
                100 * suggestive_low_ld / suggestive_snps))
message(sprintf("Genome-wide significant SNPs (P < %g): %d", P_GW, gw_snps))
message(sprintf("Genome-wide significant SNPs with low LD: %d (%.1f%%)", gw_low_ld,
                100 * gw_low_ld / gw_snps))
