# Data access

This repository contains the full analysis code for the manuscript, but not the underlying
individual-level genetic and neuropathology data, which cannot be shared publicly for privacy
reasons. This document explains what's excluded, which scripts need it, and how to obtain
equivalent data.

## What's included

- All analysis code (`data_prep/`, and the figure/table scripts at the repository root).
- Aggregate, non-identifiable inputs: GWAS summary statistics (a small window, chr6:32,000,000-
  34,000,000, around the HLA-II region - the only part Figures 1 and 2 need), Ensembl gene
  coordinates, and the small aggregate outputs of the `data_prep/` pipeline (clumping results,
  COJO haplotype table, LD matrices, regression summary statistics) — everything the 7 main
  figure/table scripts need to run as-is from a clone of this repository (see `README.md`).
- `raw_input_data/data_paths.R.example`: a template (placeholder paths only) for the private-path
  config every `data_prep/` script needs — see "How to obtain equivalent data" below.

## What's excluded

Any script reading data from this category says so in its header ("Reads private data not
included in this repository").

1. **Individual-level genetic and neuropathology data** for the three cohorts described in the
   manuscript (the study cohort, the 100-plus Study post-mortem subset, and the Netherlands Brain
   Bank replication cohort): raw genotypes/dosages, HLA imputation calls, phenotypes, and
   neuropathology/microglia measurements. This is what the `data_prep/` scripts whose header says
   "Reads private data not included in this repository" read from a local, gitignored
   `raw_input_data/` folder (or, for very large files, directly from their original location).
2. **The full genome-wide chr6 dosage file** (~27GB) that the study cohort's genetic data is
   extracted from — too large to redistribute regardless of privacy, and not specific to this
   analysis.
3. **The EADB-GWAS-2026 summary statistics**, aligned with PLINK — aggregate (not individual-level)
   data, but not ours to redistribute in full, and too large to be worth including beyond what's
   actually needed; only the chr6:32,000,000-34,000,000 window Figures 1 and 2 use is included
   (`data_prep/reference_data/gwas_summary_stats_eadb2026.txt`). The wider chromosome 6 export this
   was trimmed from is excluded like any other private input.

## How to obtain equivalent data

See the Data Availability statement in the manuscript
([doi.org/10.64898/2026.08.24.26361177](https://doi.org/10.64898/2026.08.24.26361177)): data are
available from the Alzheimer Genetics Hub (AGH, https://alzheimergenetics.org/) upon submission
of a research proposal.

If you have obtained access to this data and want to re-run the pipeline yourself, place your
copies under `raw_input_data/<cohort>/` matching the layout each `data_prep/` script's
header describes, then copy `raw_input_data/data_paths.R.example` to `raw_input_data/data_paths.R`
and fill in the placeholder paths (see `raw_input_data/SOURCE_DATA_MANIFEST.md` for what each one
is). Run every script from the repository root (see the top-level `README.md`), and run the
`data_prep/` scripts in dependency order.

Every private input this pipeline needs is one of the `data_paths.R` constants above — nothing
is read from anywhere else, and no script hardcodes a path of its own. Once every file is placed
under `raw_input_data/` as described (or, for the handful too large to copy, `data_paths.R`
points at wherever you keep them) and `data_paths.R` also points at your own PLINK2/PLINK/GCTA
install, the entire pipeline — every `data_prep/` script and all 7 main figure/table scripts —
runs end to end. Nothing else needs to change anywhere in the repository.

See [`SOFTWARE.md`](SOFTWARE.md) for exact PLINK2/PLINK/GCTA versions and the R environment
(`environment.yml`) this pipeline was run under.
