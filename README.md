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
├── DATA_ACCESS.md                # Which data are included, and how to obtain the rest
├── SOFTWARE.md                   # R, PLINK and GCTA versions and setup
├── environment.yml               # Conda environment (R and packages, version-pinned)
├── LICENSE                       # MIT license
├── CITATION.cff                  # Citation metadata
├── figure1_figure2.R             # Main scripts: produce every figure and table in the
├── figure3.R                     #   manuscript from the aggregate files in data_prep/
├── figure4_suppfigure1.R
├── figure5.R
├── table1.R
├── table2_supptable1.R
├── table3_supptable2.R
├── figures/                      # Figure outputs
├── tables/                       # Table outputs
├── data_prep/                    # Three cohort pipelines and external reference data
│   ├── study_cohort_genetics/    #   private data required (see DATA_ACCESS.md)
│   ├── 100plus_study_cohort/
│   ├── nbb_replication_cohort/
│   └── reference_data/           #   Ensembl coordinates; GWAS statistics for chr6:32-34 Mb
└── raw_input_data/               # Location for private inputs (gitignored)
    └── data_paths.R.example      #   template for raw_input_data/data_paths.R
```

The 7 main scripts read only the aggregate files produced by `data_prep/` and public reference
data, never private data. Most `data_prep/` scripts need individual-level data that is not
included here; each says so in its header (see `DATA_ACCESS.md`).

## Software and environment

The main scripts need R and a set of CRAN packages. The `data_prep/` scripts also call PLINK 1.9,
PLINK 2 and GCTA. See [`SOFTWARE.md`](SOFTWARE.md) for versions and installation. To create the R
environment:

```bash
conda env create -f environment.yml
conda activate hla-ii-finemapping
```

## Running

All paths are relative to the repository root, so run every script from there, for example
`Rscript figure1_figure2.R` or `Rscript data_prep/study_cohort_genetics/1_clumping_cojo.R`. To run
a script from elsewhere, uncomment and edit the `setwd("/path/to/repository")` line near its top.

Scripts that read individual-level data or call PLINK or GCTA also source
`raw_input_data/data_paths.R`. Copy `raw_input_data/data_paths.R.example` to
`raw_input_data/data_paths.R` and fill in the paths to your copy of the data and your PLINK and
GCTA installations (see `DATA_ACCESS.md`).

The 7 main scripts run directly on the files committed in `data_prep/reference_data/` and
`data_prep/*/output/`. Regenerating those files requires access to the cohort data (see
`DATA_ACCESS.md`), placed under `raw_input_data/` as each script's header describes. Run the
`data_prep/` pipelines in this order: `study_cohort_genetics/` first, then `100plus_study_cohort/`
(which uses step 5 of `study_cohort_genetics/`), then `nbb_replication_cohort/`. The NBB pipeline
is independent of the other two except for its last step, which also uses their outputs.

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

Most scripts run on a laptop. Two steps in `study_cohort_genetics/` need much more memory:

- `1_clumping_cojo.R` runs PLINK 2 extraction and clumping on the ~27 GB genome-wide
  `CHR6_GENOMEWIDE_DOSAGE` file and exports a per-sample dosage file, `chr6_mhc.raw.gz`
  (~5.7 GB before compression, ~360 MB after).
- `5_get_snps_alleles_dosages.R` reads that file back in with `fread()`.

Both need considerably more memory than the size of their input. On a shared HPC system, submit
them as batch jobs rather than running them on a login node.

## Archival

A permanent copy of this code is archived on Zenodo: [doi:10.5281/zenodo.23219197](https://doi.org/10.5281/zenodo.23219197). This DOI always resolves to the latest release.

## License

Released under the [MIT License](LICENSE).

## Contact

For questions about the code or analysis, contact the corresponding author:
Henne Holstege, h.holstege@amsterdamumc.nl
