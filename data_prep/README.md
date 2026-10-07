# data_prep/

Three cohort pipelines, each turning private individual-level data into the aggregate files that
the 7 main figure and table scripts at the repository root use, plus a folder of external
reference data:

- `study_cohort_genetics/`: the genetic cohort of ~6,053 individuals (CHCs, controls and AD)
  underlying Figures 1-4.
- `100plus_study_cohort/`: the 100-plus Study post-mortem subset (89 CHCs, 7 AD), for Figure 5A/B
  and Table 3.
- `nbb_replication_cohort/`: the independent Netherlands Brain Bank replication cohort, for
  Figure 5C and Table 3.
- `reference_data/`: external reference files with no generating script. Ensembl gene
  coordinates, and the chr6:32,000,000-34,000,000 window of the EADB-GWAS-2026 summary statistics
  used for Figures 1 and 2 (see `DATA_ACCESS.md`).

Scripts are numbered in the order they run within each pipeline, and each header documents its
inputs and outputs. Scripts whose header says "Reads private data not included in this
repository" need individual-level data (see `DATA_ACCESS.md`).

Each script writes to its pipeline's `output/` folder. Most of these files are individual-level
data or large intermediates (PLINK and GCTA working files, genotype exports) and are not
committed; only a few small aggregate files are (see `.gitignore`). The combined SNP, allele and
neuropathology files of the 100-plus and NBB cohorts are individual-level, so a README in each of
those `output/` folders explains why they are missing.
