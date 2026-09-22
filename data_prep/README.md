# data_prep/

Three cohort-specific pipelines, each turning private individual-level data into the small
aggregate files the 7 main figure/table scripts at the repository root need, plus one folder of
static external reference data:

- `study_cohort_genetics/` — the ~6,053-person genetic cohort (CHCs, controls, AD) underlying
  Figures 1-4.
- `100plus_study_cohort/` — the 100-plus Study post-mortem subset (89 CHCs, 7 AD), for Figure 5A/B
  and Table 3.
- `nbb_replication_cohort/` — the independent Netherlands Brain Bank replication cohort, for
  Figure 5C and Table 3.
- `reference_data/` — external reference files with no generating script: Ensembl gene
  coordinates (committed as-is), and a small chr6:32,000,000-34,000,000 window of the
  EADB-GWAS-2026 summary statistics, trimmed from a much larger private export to the region
  Figures 1 and 2 actually need - see `DATA_ACCESS.md`.

Each script in the three cohort pipelines has a numbered filename giving its position in that
cohort's pipeline, and a header documenting its inputs/outputs. Scripts whose header says "Reads
private data not included in this repository" read individual-level data not included here — see
`DATA_ACCESS.md` at the repository root.

Every cohort script writes into its own `output/` subfolder. Most of what lands there is either raw
individual-level data or large regenerable intermediates (PLINK/GCTA working files, full genotype
exports), so it isn't committed — only a handful of small, aggregate, shareable files are (see the
repository's `.gitignore`). That's why `study_cohort_genetics/output/` shows only a few files despite
the scripts writing many more, and why `100plus_study_cohort/output/` and
`nbb_replication_cohort/output/`'s combined SNP+allele+neuropathology files contain nothing but a
placeholder `README.md` explaining the same thing — those two are individual-level.
