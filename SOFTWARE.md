# Software requirements

This document covers all software needed to run this pipeline: the `data_prep/` scripts
(private data required — see `DATA_ACCESS.md`) and the 7 main figure/table scripts at the
repository root (run directly against the small aggregate files already committed here).

---

## R environment

Every script was run under the conda environment defined in `environment.yml`. It pins
R 4.6.1 and the exact version of every CRAN package any script in this repository loads.

### Recreating the environment

```bash
conda env create -f environment.yml
conda activate hla-ii-finemapping
```

### Which scripts need which packages

The 7 main scripts at the repository root, and `data_prep/`, together load: `data.table`,
`dplyr`, `stringr`, `tidyr`, `ggplot2`, `patchwork`, `RColorBrewer`, `gridExtra`, `ggpubr`,
`corrplot`, `readxl` and `scales` (the last used inline as `scales::alpha()`, never
`library()`-loaded). No single script needs all of them — see each script's own header for
its specific inputs, and its `library()` calls for what it actually loads.

---

## PLINK

Two versions of PLINK are used at different pipeline steps, both configured in
`raw_input_data/data_paths.R` (`PLINK1_PATH`, `PLINK2_PATH`).

### PLINK 2 (`plink2`)

Used for almost everything: region extraction (`--make-pgen`), clumping (`--clump`),
per-sample dosage export (`--export A`), format conversion (`--make-bed`), and the
`--import-dosage` step that builds Figure 3's LD matrix.

Download: [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/)

Version used in this study: **v2.00a6LM AVX2 Intel (18 Aug 2024)**.

### PLINK 1.9 (`plink`)

Used only for `--r2` (Figure 1's regional LD) and `--r2 square` (Figure 3's LD matrix) —
PLINK 2 doesn't support `--r2` the same way.

Download: [https://www.cog-genomics.org/plink/1.9/](https://www.cog-genomics.org/plink/1.9/)

Version used in this study: **v1.90b6.26 (2 Apr 2022)**.

---

## GCTA

Used by `data_prep/study_cohort_genetics/1_clumping_cojo.R` for the Conditional & Joint
(COJO) stepwise analysis (`--cojo-slct`) that identifies the independent haplotype lead
SNPs reported in Table 1, configured via `GCTA_PATH` in `raw_input_data/data_paths.R`.

- Publication: Yang et al. (2011) *American Journal of Human Genetics*
  [doi:10.1016/j.ajhg.2010.11.011](https://doi.org/10.1016/j.ajhg.2010.11.011)
- Download: [https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)

Version used in this study: **v1.94.1**.

---

## HLA imputation

`raw_input_data/*/hla_imputation/result_{DQA1,DQB1,DRB1}.txt` (both cohorts) are HIBAG output
— classical two-field HLA allele calls imputed from SNP array genotypes. Generating them is
necessary before this pipeline can run, but imputation itself is not addressed by any script
in this repository: the results are supplied as private raw input (see `DATA_ACCESS.md` and
`raw_input_data/SOURCE_DATA_MANIFEST.md`).

For the imputation pipeline itself (HIBAG model, liftover, PLINK steps), see
[`hla_imputation/`](https://github.com/dalvarez15/Carrying-specific-HLA-alleles-decreases-the-chance-of-reaching-healthy-old-age/tree/main/hla_imputation)
in the code repository for our earlier publication (Álvarez Sirvent et al., *npj Aging*, 2026),
which uses the same imputation approach.
