# Integration Tests for prolfquapp CLI Pipelines

Cross-package integration tests for `prolfqua_dea.sh` (CMD_DEA_V2.R) and `prolfqua_qc.sh` (CMD_QUANT_QC.R). These test the full pipeline end-to-end using small fixture datasets (~100 proteins) subsetted from real data.

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

## Repository and Data Policy

The Git repository should stay small enough to clone and review quickly. Track source code, test code, scripts,
documentation, small configs, and fixture metadata. Do not track generated outputs or bulky real-data payloads.

Tracked in Git:
- `README.md`, `CLAUDE.md`, `Makefile`, `.gitignore`
- `R/`, `scripts/`, `tests/`, and TODO/planning documents
- lightweight fixture documentation such as `fixtures/README.md`

Ignored locally:
- `fixtures/*` payload directories
- `.cache/`
- `logs/`
- `test-outputs/`
- `reference/`
- `prolfquapp_docker.sh`

The current local checkout can be multiple gigabytes because `test-outputs/` duplicates inputs and rendered reports.
This is expected locally, but those files should not be committed. If a large fixture becomes required for reproducible
remote use, prefer an explicit download/regeneration step with checksums. Use Git LFS only after deciding that a fixture
must live with the repository despite its size.

When publishing this project, a natural remote is a dedicated repository under the `prolfqua` GitHub organization. Before
the first push, check the candidate commit with:

```bash
git status --short --ignored
git ls-files -z | xargs -0 du -h | sort -h | tail
```

Only ignored local data should account for large disk usage.

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
make test-dea-internal
# etc. — see make help for the full list
```

### Run the WU345302 facade matrix

The WU345302 facade matrix is a broader CLI smoke test for all registered DEA facade models on the same DIA-NN fixture.
It is useful after changing model facades, contrast adapters, or peptide-to-protein workflows.

```bash
make wu345302-facades
```

This target runs:

```bash
bash scripts/run_wu345302_facades.sh
```

The script writes its outputs under `test-outputs/wu345302_facades/`:

| File | Purpose |
|------|---------|
| `status.tsv` | Per-model command status and exit code |
| `failures.tsv` | Failed models with log file paths and extracted error messages |
| `model_summary.tsv` | Row counts, finite FDR/diff counts, significant counts, and report links |
| `pairwise_vs_limma_impute.tsv` | Pairwise comparisons against the `limma_impute` reference model |

The facade matrix uses installed package entry points by default. Run `make install` first when you want the installed
CLI scripts to reflect the current local package sources.

## What's in here

```
integration_test/
  R/
    create_test_fixtures.R         # One-time fixture generator
  fixtures/
    README.md                      # Fixture data policy; payload directories are ignored
  reference/                       # Ignored regression references after make save-references
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
      test-dea-internal-calibration.R # DEA with internal calibration
  scripts/
    run_wu345302_facades.sh        # Run all registered WU345302 facade models
    summarize_wu345302_facades.R   # Summarize model outputs and write failures.tsv
```

## Test Details

| Test file | Script tested | Software flag | Fixture | What it checks |
|-----------|--------------|---------------|---------|----------------|
| test-dea-maxquant | CMD_DEA_V2.R | `prolfquapp.MAXQUANT` | maxquant_ionstar | HTML, XLSX, SE.rds, RNK, parquet; SE has contrast columns with diff/FDR |
| test-dea-msstats | CMD_DEA_V2.R | `prolfquapp.MSSTATS` | fragpipe_ionstar | Same outputs; different preprocessor path |
| test-dea-fp-tmt | CMD_DEA_V2.R | `prolfquapp.FP_TMT` | fp_tmt_total | Same + verifies >=4 complex contrasts (2x3 factorial design) |
| test-dea-fp-singlesite | CMD_DEA_V2.R | `prolfquappPTMreaders.FP_singlesite` | fp_singlesite_phospho | Same + PTM site-level aggregation; skips if prolfquappPTMreaders not installed |
| test-qc-maxquant | CMD_QUANT_QC.R | `MAXQUANT` | maxquant_ionstar | HTML reports + XLSX produced |
| test-dea-regression | CMD_DEA_V2.R | all | all fixtures | Compares SE fold-changes against saved references (correlation > 0.999) |
| test-dea-internal | CMD_DEA_V2.R | internal fixture config | internal calibration fixture | Checks `center_to_reference` internal calibration workflow |

### Software naming gotcha

CMD_DEA_V2.R uses `get_procfuncs()` which returns **prefixed** keys like `prolfquapp.MAXQUANT`. CMD_QUANT_QC.R uses `prolfqua_preprocess_functions` directly with **unprefixed** keys like `MAXQUANT`. The test files and fixture configs reflect this difference.

## Regression Tests

The regression test (`test-dea-regression.R`) compares SummarizedExperiment outputs against saved baselines. **References must be generated from the released Docker image** (`prolfqua/prolfquapp:latest`) to detect regressions against the published version.

```bash
make save-references-docker   # generate baselines from released Docker image
make test-dea-regression      # compare local dev output against Docker baselines
```

On ARM Mac, Docker runs with `--platform linux/amd64` (Rosetta emulation). The `fp_singlesite_phospho` fixture is skipped because the Docker image doesn't include `prolfquappPTMreaders`.

A local `save-references` target also exists for re-baselining after a new release:

```bash
make save-references    # generate baselines from local dev install
```

## How it works

Each test:

1. Copies a fixture directory into a fresh temp dir
2. Runs `CMD_DEA_V2.R` or `CMD_QUANT_QC.R` via `system()` (as a subprocess, same as the shell scripts)
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
| diann_wu345302 | local fixture payload | DIA-NN report, FASTA files, dataset, and config template | WU345302 facade matrix |

## Regenerating fixtures

```bash
make clean              # remove generated fixtures and outputs, keep diann_wu345302
make fixtures           # regenerate from source data
```

`make fixtures` regenerates the standard MaxQuant, MSstats, FP_TMT, and FP_singlesite fixtures from `prolfquadata` and
the Zenodo PTM archive. It does not currently recreate the local `diann_wu345302` payload used by
`make wu345302-facades`.

`make clean` removes the generated standard fixtures, logs, and development outputs, but intentionally keeps
`fixtures/diann_wu345302` because that local DIA-NN payload is not recreated by `make fixtures`.
