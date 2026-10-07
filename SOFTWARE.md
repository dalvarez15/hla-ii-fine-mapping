# Software requirements

This document lists the software needed to run the `data_prep/` scripts (private data required,
see `DATA_ACCESS.md`) and the 7 main figure and table scripts at the repository root.

---

## R environment

The analysis was run with R 4.6.1. `environment.yml` provides R 4.5.3 with the same package
versions apart from minor patch releases (data.table 1.18.4, gridExtra 2.3.1, readxl 1.5.0.1); it
reproduces the figures and tables in this repository. Figure 5 is saved with `cairo_pdf()`, so its
font depends on the fonts installed on the system.

### Creating the environment

```bash
conda env create -f environment.yml
conda activate hla-ii-finemapping
```

### Packages

The scripts use `data.table`, `dplyr`, `stringr`, `tidyr`, `ggplot2`, `patchwork`,
`RColorBrewer`, `gridExtra`, `ggpubr`, `corrplot`, `readxl` and `scales`. Each script loads only
the packages it needs.

---

## PLINK

Two versions of PLINK are used, configured in `raw_input_data/data_paths.R` as `PLINK2_PATH` and
`PLINK1_PATH`.

### PLINK 2 (`plink2`)

Used for region extraction (`--make-pgen`), clumping (`--clump`), per-sample dosage export
(`--export A`), format conversion (`--make-bed`) and the `--import-dosage` step that builds the
LD matrix for Figure 3.

Download: [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/)

Version used in this study: **v2.00a6LM AVX2 Intel (18 Aug 2024)**.

### PLINK 1.9 (`plink`)

Used only for `--r2` (regional LD, Figure 1) and `--r2 square` (LD matrix, Figure 3), which
PLINK 2 does not provide in the same form.

Download: [https://www.cog-genomics.org/plink/1.9/](https://www.cog-genomics.org/plink/1.9/)

Version used in this study: **v1.90b6.26 (2 Apr 2022)**.

---

## GCTA

Used by `data_prep/study_cohort_genetics/1_clumping_cojo.R` for the conditional and joint (COJO)
stepwise analysis (`--cojo-slct`) that identifies the independent haplotype lead SNPs in Table 1.
Configured as `GCTA_PATH` in `raw_input_data/data_paths.R`.

- Publication: Yang et al. (2011) *American Journal of Human Genetics*
  [doi:10.1016/j.ajhg.2010.11.011](https://doi.org/10.1016/j.ajhg.2010.11.011)
- Download: [https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)

Version used in this study: **v1.94.1**.

---

## HLA imputation

The files `raw_input_data/*/hla_imputation/result_{DQA1,DQB1,DRB1}.txt` (both cohorts) are HIBAG
output: classical two-field HLA allele calls imputed from SNP array genotypes. They are required
inputs, but the imputation itself is not part of this repository; the results are supplied as
private raw input (see `DATA_ACCESS.md`).

The imputation pipeline (HIBAG model, liftover, PLINK steps) is in
[`hla_imputation/`](https://github.com/dalvarez15/Carrying-specific-HLA-alleles-decreases-the-chance-of-reaching-healthy-old-age/tree/main/hla_imputation)
in the repository for our earlier publication (Álvarez Sirvent et al., *npj Aging*, 2026), which
uses the same imputation approach.
