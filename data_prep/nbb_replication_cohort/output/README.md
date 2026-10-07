This folder is filled by running `../1_snps_phenotypes.R` and then
`../2_combine_snps_alleles_neuropathology_nbb.R` with access to the individual-level data
described in their headers. Their output, `hla_alleles_snps_dosages_neuropathology_nbb.csv`, has
one row per participant, so it is not committed to this repository (see `DATA_ACCESS.md` at the
repository root).

`microglia_neuropathology_regression.csv`, written by `../3_microglia_neuropathology_regression.R`,
contains aggregate association statistics and is committed.
