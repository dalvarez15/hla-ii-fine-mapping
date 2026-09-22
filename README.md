# HLA-II haplotype fine-mapping

Code repository for:

> Álvarez Sirvent D, Luimes MC, Tesi N, Rohde SK, Salazar AN, Bijker LJ, van Schoor NM, Tijms BM,
> Vijverberg EGB, Strijbis EMM, Hoekstra EJ, Holtman IR, van der Lee SJ, Hulsman M, Holstege H.
> *Fine-mapping HLA-II haplotypes in Alzheimer's disease and healthy longevity reveals distinct
> associations with microglial HLA-II load and neuropathology.* medRxiv 2026.08.24.26361177
> (2026). [doi.org/10.64898/2026.08.24.26361177](https://doi.org/10.64898/2026.08.24.26361177)

---

## Repository structure

```
.
├── README.md                     # This file
├── DATA_ACCESS.md                # What's public vs. private, and how to obtain the private data
├── SOFTWARE.md                   # R, PLINK, GCTA versions and setup
├── environment.yml               # Conda environment (R + all packages, version-pinned)
├── figure1_figure2.R             # Main scripts: produce every figure/table in the manuscript
├── figure3.R                     #   from the small aggregate files data_prep/ produces
├── figure4_suppfigure1.R         #
├── figure5.R                     #
├── table1.R                      #
├── table2_supptable1.R           #
├── table3_supptable2.R           #
├── figures/                      # Script outputs (figures)
├── tables/                       # Script outputs (tables)
├── data_prep/                    # Three cohort-specific pipelines + external reference data
│   ├── study_cohort_genetics/    #   private data required - see DATA_ACCESS.md
│   ├── 100plus_study_cohort/     #
│   ├── nbb_replication_cohort/   #
│   └── reference_data/           #   Ensembl coords as-is; GWAS stats trimmed to a small window
└── raw_input_data/               # Private inputs go here (gitignored) - see DATA_ACCESS.md
    └── data_paths.R.example      #   template for raw_input_data/data_paths.R
```

The 7 root scripts read only the small aggregate files `data_prep/` produces (or reference data
with no privacy restrictions) — none of them read private data directly. Most `data_prep/` scripts
require individual-level data not included in this repository; each says so in its own header —
see `DATA_ACCESS.md`.

## Software and environment

This pipeline needs R plus a handful of CRAN packages (all 7 main scripts), and additionally
PLINK 1.9, PLINK 2 and GCTA for the `data_prep/` scripts that call them. See
[`SOFTWARE.md`](SOFTWARE.md) for exact versions and install instructions. To recreate the R
environment used for this analysis:

```bash
conda env create -f environment.yml
conda activate hla-ii-finemapping
```

## Running

Every script uses paths relative to the repository root, so run them **from the repository root**,
e.g. `Rscript figure1_figure2.R` or `Rscript data_prep/study_cohort_genetics/1_clumping_cojo.R` — not
from inside `data_prep/` or with a different working directory. Alternatively, uncomment and edit
the `setwd("/path/to/repository")` line near the top of the script to set the working directory
explicitly, which lets you run it from anywhere. Scripts that read individual-level data or call
external tools (PLINK, GCTA) also source `raw_input_data/data_paths.R`; copy
`raw_input_data/data_paths.R.example` to `raw_input_data/data_paths.R` and fill in the paths for
your own copy of the data and your own PLINK/GCTA install — see `DATA_ACCESS.md` and
`raw_input_data/SOURCE_DATA_MANIFEST.md` for what each path is.

All 7 main scripts can be run directly against the files already committed in this repository
(`data_prep/reference_data/` and `data_prep/*/output/`). To regenerate those files from scratch,
you need access to the underlying cohort data — see `DATA_ACCESS.md` — placed under
`raw_input_data/` as each script's header describes, then run the `data_prep/` scripts in order:
`study_cohort_genetics/` first, then `100plus_study_cohort/` (which needs `study_cohort_genetics/`'s
step 5), then `nbb_replication_cohort/` (independent of the other two cohorts, except for its last
step, which combines the association statistics from all three):

```
data_prep/study_cohort_genetics/1_clumping_cojo.R
data_prep/study_cohort_genetics/2_regional_ld.R
data_prep/study_cohort_genetics/3_cojo_annotate.R
data_prep/study_cohort_genetics/4_read_hla_imputation.R
data_prep/study_cohort_genetics/5_get_snps_alleles_dosages.R
data_prep/study_cohort_genetics/6_haplotype_ld.R
data_prep/study_cohort_genetics/7_allele_snp_regression.R

data_prep/100plus_study_cohort/1_combine_snps_alleles_neuropathology_100plus.R

data_prep/nbb_replication_cohort/1_snps_phenotypes.R
data_prep/nbb_replication_cohort/2_combine_snps_alleles_neuropathology_nbb.R
data_prep/nbb_replication_cohort/3_microglia_neuropathology_regression.R
```

### Resource requirements

Most scripts are lightweight and run fine on a laptop. Two steps in `study_cohort_genetics/` are
not:

- `1_clumping_cojo.R` extracts from and runs PLINK2 clumping against the ~27GB genome-wide
  `CHR6_GENOMEWIDE_DOSAGE` file, and exports a per-sample dosage file (`chr6_mhc.raw.gz`,
  ~5.7GB before it gzips it in place - about 360MB compressed).
- `5_get_snps_alleles_dosages.R` reads that file back in with `fread()`.

Both need substantially more memory than the input file size to run comfortably, and on a shared
HPC system should be submitted as a batch job rather than run on a login node, which is typically
memory-limited and shared across many users. Every other script operates on much smaller,
already-trimmed data and doesn't need special handling.

## Contact

For questions about the code or analysis, contact the corresponding author:
Henne Holstege — h.holstege@amsterdamumc.nl
