# Integration Tests for prolfquapp CLI Pipelines

Cross-package integration tests for `prolfqua_dea.sh` (CMD_DEA.R) and `prolfqua_qc.sh` (CMD_QUANT_QC.R). These test the full pipeline end-to-end using small fixture datasets (~100 proteins) subsetted from real data.

## Prerequisites

The fixture generator needs:
- `prolfquadata` package installed (IonStar MaxQuant + MSFragger ZIPs)
- `prophosqua` repo checked out alongside the other packages (PTM example data)
- `seqinr` package

## Quick Start

All commands run from `integration_test/`:

```bash
cd integration_test
make help               # show all targets
make install            # reinstall prolfqua + prolfquapp + prolfquappPTMreaders from local source
make fixtures           # generate fixture data from real datasets (one-time, ~3 min)
make test               # run all integration tests (~2-3 min)
```

### Typical workflow after editing prolfqua

```bash
cd integration_test
make install            # reinstall packages with your changes
make test               # check nothing broke
```

### Run a single test

```bash
make test-dea-maxquant
make test-qc-maxquant
make test-dea-fp-singlesite
# etc. — see make help for the full list
```

## What's in here

```
integration_test/
  R/
    create_test_fixtures.R         # One-time fixture generator
  fixtures/                        # Generated test data (not checked in)
    maxquant_ionstar/              # ~100 proteins, peptides.txt format
    fragpipe_ionstar/              # ~50 proteins, MSstats.csv format
    fp_tmt_total/                  # ~100 proteins, psm.tsv (TMT)
    fp_singlesite_phospho/         # ~50 proteins, abundance_single-site_None.tsv
    reference/                     # Baseline SummarizedExperiment.rds (after make save-references)
  tests/
    testthat.R                     # Entry point for testthat::test_dir()
    testthat/
      helper-common.R              # Shared utilities (run_dea, run_qc, find_*_outputs)
      test-dea-maxquant.R          # DEA with MAXQUANT preprocessor
      test-dea-msstats.R           # DEA with MSSTATS preprocessor
      test-dea-fp-tmt.R            # DEA with FP_TMT preprocessor (TMT, VSN, complex contrasts)
      test-dea-fp-singlesite.R     # DEA with FP_singlesite preprocessor (phospho PTM)
      test-qc-maxquant.R           # QC pipeline (CMD_QUANT_QC.R)
      test-dea-regression.R        # Compare outputs against saved reference SE objects
```

## Test Details

| Test file | Script tested | Software flag | Fixture | What it checks |
|-----------|--------------|---------------|---------|----------------|
| test-dea-maxquant | CMD_DEA.R | `prolfquapp.MAXQUANT` | maxquant_ionstar | HTML, XLSX, SE.rds, RNK, parquet; SE has contrast columns with diff/FDR |
| test-dea-msstats | CMD_DEA.R | `prolfquapp.MSSTATS` | fragpipe_ionstar | Same outputs; different preprocessor path |
| test-dea-fp-tmt | CMD_DEA.R | `prolfquapp.FP_TMT` | fp_tmt_total | Same + verifies >=4 complex contrasts (2x3 factorial design) |
| test-dea-fp-singlesite | CMD_DEA.R | `prolfquappPTMreaders.FP_singlesite` | fp_singlesite_phospho | Same + PTM site-level aggregation; skips if prolfquappPTMreaders not installed |
| test-qc-maxquant | CMD_QUANT_QC.R | `MAXQUANT` | maxquant_ionstar | HTML reports + XLSX produced |
| test-dea-regression | CMD_DEA.R | all | all fixtures | Compares SE fold-changes against saved references (tolerance 1e-3) |

### Software naming gotcha

CMD_DEA.R uses `get_procfuncs()` which returns **prefixed** keys like `prolfquapp.MAXQUANT`. CMD_QUANT_QC.R uses `prolfqua_preprocess_functions` directly with **unprefixed** keys like `MAXQUANT`. The test files and fixture configs reflect this difference.

## Regression Tests

The regression test (`test-dea-regression.R`) compares SummarizedExperiment outputs against saved baselines. References don't exist initially — tests skip with a message.

```bash
make save-references    # run all 4 DEA pipelines, save baseline SE objects
make test-dea-regression  # compare new runs against baselines
```

After intentional changes to the analysis results:

```bash
make clean-references   # remove old baselines
make save-references    # regenerate
```

## How it works

Each test:

1. Copies a fixture directory into a fresh temp dir
2. Runs `CMD_DEA.R` or `CMD_QUANT_QC.R` via `system()` (as a subprocess, same as the shell scripts)
3. Asserts expected output files exist
4. For DEA: loads `SummarizedExperiment.rds` and checks structure (dimensions, nested contrast DataFrames with `diff`/`FDR` columns, non-NA values)
5. Cleans up the temp dir

The helper functions `run_dea()` and `run_qc()` in `helper-common.R` handle the temp dir setup, argument construction, and subprocess invocation. They use `cd workdir && Rscript ...` so that any files the scripts write to the current directory land in the temp dir.

## Fixture Data Sources

| Fixture | Source package | Original file | Subset size |
|---------|---------------|---------------|-------------|
| maxquant_ionstar | prolfquadata | IonStar MaxQuant ZIP / `peptides.txt` | ~100 proteins, groups b vs e |
| fragpipe_ionstar | prolfquadata | IonStar MSFragger ZIP / `MSstats.csv` | ~50 proteins, groups B vs E |
| fp_tmt_total | prophosqua | PTM_example_analysis_v2 / `psm.tsv` (70 MB) | ~100 proteins, 22 samples, 4 contrasts |
| fp_singlesite_phospho | prophosqua | PTM_example_analysis_v2 / `abundance_single-site_None.tsv` | ~50 proteins (multiple sites each), 22 samples, 4 contrasts |

## Regenerating fixtures

```bash
make clean              # remove all fixtures
make fixtures           # regenerate from source data
```
