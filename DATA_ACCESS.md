# Data access

This repository contains the full analysis code for the manuscript, but not the individual-level
genetic and neuropathology data, which cannot be shared publicly for privacy reasons. This
document lists what is excluded, which scripts need it, and how to obtain it.

## What's included

- All analysis code (`data_prep/` and the figure and table scripts at the repository root).
- Aggregate, non-identifiable inputs: GWAS summary statistics for chr6:32,000,000-34,000,000
  around the HLA-II region (the window Figures 1 and 2 use), Ensembl gene coordinates, and the
  aggregate outputs of the `data_prep/` pipelines (clumping results, COJO haplotype table, LD
  matrices and regression summary statistics). These are all the 7 main figure and table scripts
  need to run from a clone of this repository (see `README.md`).
- `raw_input_data/data_paths.R.example`: a template, with placeholder paths only, for the path
  configuration the `data_prep/` scripts need.

## What's excluded

Scripts that read these data say so in their header ("Reads private data not included in this
repository").

1. **Individual-level genetic and neuropathology data** for the three cohorts in the manuscript
   (the study cohort, the 100-plus Study post-mortem subset and the Netherlands Brain Bank
   replication cohort): genotypes and dosages, HLA imputation calls, phenotypes, and
   neuropathology and microglia measurements. The `data_prep/` scripts read these from a local,
   gitignored `raw_input_data/` folder, or for very large files from their original location.
2. **The genome-wide chr6 dosage file** (~27 GB) from which the study cohort's genetic data are
   extracted. It is too large to redistribute and not specific to this analysis.
3. **The EADB-GWAS-2026 summary statistics**, aligned with PLINK. These are aggregate data, but
   not ours to redistribute in full. Only the chr6:32,000,000-34,000,000 window used for Figures 1
   and 2 is included (`data_prep/reference_data/gwas_summary_stats_eadb2026.txt`).

## How to obtain the data

See the Data Availability statement in the manuscript
([doi.org/10.64898/2026.08.24.26361177](https://doi.org/10.64898/2026.08.24.26361177)): data are
available from the Alzheimer Genetics Hub (AGH, https://alzheimergenetics.org/) upon submission
of a research proposal.

With access to the data, place your copies under `raw_input_data/<cohort>/` following the layout
in each `data_prep/` script's header. Copy `raw_input_data/data_paths.R.example` to
`raw_input_data/data_paths.R` and fill in the placeholder paths; the comments in the template
describe each one. Files too large to copy can stay where they are, with `data_paths.R` pointing
to them. The PLINK 2, PLINK 1.9 and GCTA paths are set in the same file. All private inputs are
read through `data_paths.R`, so nothing else in the repository needs to change. Run every script
from the repository root, and the `data_prep/` scripts in the order given in `README.md`.

See [`SOFTWARE.md`](SOFTWARE.md) for tool versions and the R environment.
