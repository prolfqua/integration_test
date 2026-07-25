# Integration Tests for prolfquapp CLI Pipelines

Cross-package integration tests for `prolfqua_dea.sh` (CMD_DEA_V2.R) and `prolfqua_qc.sh` (CMD_QUANT_QC.R). These test
the full pipeline end-to-end using small fixture datasets derived from real data.

## Prerequisites

The fixture generator needs:
- `prolfquadata` package installed (IonStar MaxQuant + MSFragger ZIPs)
- `prolfquappPTMreaders` installed (bundled Spectronaut BGS example)
- `seqinr` package
- an Internet connection for the cached FragPipe PTM example download

## Quick Start

All commands run from `integration_test/`:

```bash
cd integration_test
make help               # show all targets
make install            # reinstall prolfqua + prolfquapp + prolfquappPTMreaders from local source
make fixtures           # generate fixture data from real datasets (one-time, ~3 min)
make test               # run all integration tests (~10-12 min)
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
make test-dea-fp-multisite
make test-dea-fp-combined-sty
make test-dea-bgs-site
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
  tests/
    testthat.R                     # Entry point for testthat::test_dir()
    testthat/
      helper-common.R              # Shared utilities (run_dea, run_qc, find_*_outputs)
      test-dea-maxquant.R          # DEA with MAXQUANT preprocessor
      test-dea-msstats.R           # DEA with MSSTATS preprocessor
      test-dea-fp-tmt.R            # DEA with FP_TMT preprocessor (TMT, VSN, complex contrasts)
      test-dea-fp-singlesite.R     # DEA with FP_singlesite preprocessor (phospho PTM)
      test-dea-fp-multisite.R      # DEA with FP_multisite preprocessor (phospho PTM)
      test-dea-fp-combined-sty.R   # DEA with FP_combined_STY preprocessor (phospho PTM)
      test-dea-bgs-site.R          # DEA with Spectronaut BGS_site preprocessor (phospho PTM)
      test-qc-maxquant.R           # QC pipeline (CMD_QUANT_QC.R)
      test-dea-internal-calibration.R # DEA with internal calibration
  scripts/
    docker_regression.sh           # Docker model-matrix regression driver
    compare_docker_versions.R      # Correlate DEA outputs across runs (the comparator)
    run_wu345302_facades.sh        # Run all registered WU345302 facade models
    summarize_wu345302_facades.R   # Summarize model outputs and write failures.tsv
```

## Test Details

| Test file | Script tested | Software flag | Fixture | What it checks |
|-----------|--------------|---------------|---------|----------------|
| test-dea-maxquant | CMD_DEA_V2.R | `prolfquapp.MAXQUANT` | maxquant_ionstar | HTML, XLSX, SE.rds, RNK, parquet; SE has contrast columns with diff/FDR |
| test-dea-msstats | CMD_DEA_V2.R | `prolfquapp.MSSTATS` | fragpipe_ionstar | Same outputs; different preprocessor path |
| test-dea-fp-tmt | CMD_DEA_V2.R | `prolfquapp.FP_TMT` | fp_tmt_total | Same + verifies >=4 complex contrasts (2x3 factorial design) |
| test-dea-fp-singlesite | CMD_DEA_V2.R | `prolfquappPTMreaders.FP_singlesite` | fp_singlesite_phospho | Same + PTM site-level aggregation and complete protein annotation |
| test-dea-fp-multisite | CMD_DEA_V2.R | `prolfquappPTMreaders.FP_multisite` | fp_multisite_phospho | Same + multisite hierarchy and complete protein annotation |
| test-dea-fp-combined-sty | CMD_DEA_V2.R | `prolfquappPTMreaders.FP_combined_STY` | fp_combined_sty_phospho | Same + combined-STY hierarchy and complete protein annotation |
| test-dea-bgs-site | CMD_DEA_V2.R | `prolfquappPTMreaders.BGS_site` | bgs_site_phospho | Same + Spectronaut site hierarchy and complete protein annotation |
| test-qc-maxquant | CMD_QUANT_QC.R | `MAXQUANT` | maxquant_ionstar | HTML reports + XLSX produced |
| test-dea-internal | CMD_DEA_V2.R | internal fixture config | internal calibration fixture | Checks `center_to_reference` internal calibration workflow |

These functional/feature tests run the **locally installed** prolfquapp (run `make install` first). Numeric
**regression** is handled separately and is **Docker-only** — see below.

### Software naming gotcha

CMD_DEA_V2.R uses `get_procfuncs()` which returns **prefixed** keys like `prolfquapp.MAXQUANT`. CMD_QUANT_QC.R uses `prolfqua_preprocess_functions` directly with **unprefixed** keys like `MAXQUANT`. The test files and fixture configs reflect this difference.

## Docker regression

Numeric regression is **Docker-only**: it runs a fixture through *released* `prolfqua/prolfquapp`
images (never the local dev install) and correlates the protein-level fold-change (`diff`), `p.value`
and `FDR`. There are no committed reference `.rds` blobs and no local-dev baselines — each run is
reproduced on the fly, so a check can never drift against a stale binary.

```bash
make regression        # maxquant_ionstar model matrix; fails if any r < 0.95
```

`make regression` calls `scripts/docker_regression.sh`, which runs the model matrix and then
`scripts/compare_docker_versions.R` (reads each `DE_*.xlsx` `diff_exp_analysis` sheet + the
`diff_exp_analysis_wide` FDR columns, joins on `protein_Id` × `contrast`, reports Pearson `r` and max
absolute difference, exits non-zero below `--min-cor`).

### Reference chain and model matrix

| Tag | Role | Models | How the model is chosen |
|-----|------|--------|-------------------------|
| `0.1.8` | old reference | `lm`, `lm_missing` only | config `model: prolfqua` + `model_missing:` (no `-m`) |
| `2.2.6` | model-complete reference | full facade set | `-m <facade>` overrides the config |

The matrix correlates the models that overlap, plus the new rescue against the legacy one:

| Comparison | What it guards |
|------------|----------------|
| `0.1.8/lm` ↔ `2.2.6/lm` | the core linear model across the 14-month span |
| `0.1.8/lm_missing` ↔ `2.2.6/lm_missing` | the legacy missing-group model |
| `2.2.6/lm_impute` ↔ `2.2.6/lm_missing` | the new LOD-imputation rescue vs the legacy approach |

To run another fixture / version chain, call the driver directly, e.g.
`bash scripts/docker_regression.sh fragpipe_ionstar prolfquapp.MSSTATS dataset.csv 0.95`.

### CLI differences across generations

The CLI is **not** stable across the full version range, so the runs cannot share one config:

| Generation | DEA entrypoint | `-s` software key | Config schema |
|------------|----------------|-------------------|---------------|
| `0.1.x`    | `CMD_DEA.R`     | unprefixed (`MAXQUANT`, `DIANN`, `FP_TMT`, `MSSTATS`) | older R6 config — **rejects** the 2.x YAML (`cannot add bindings to a locked environment`) |
| `2.0.x` / `2.2.x` | `CMD_DEA_V2.R` | prefixed (`prolfquapp.MAXQUANT`) | current YAML schema |

For `0.1.x` the config must be generated with that image's own `prolfqua_yaml.sh`. Conveniently its
defaults (`model=prolfqua`, `aggregate=medpolish`, `transform=robscale`, `FDR=0.1`) match the 2.x
fixture config, so the comparison is not confounded by different processing options.

### Manual cross-version comparison

`make regression` covers the model matrix on one fixture. To extend the chain to more tags (the default
`model: prolfqua` resolution per version) by hand:

```bash
# 0.1.x: generate a native config, then run with the unprefixed software key
mkdir -p test-outputs/vcompare/v018 && cp -r fixtures/maxquant_ionstar/* test-outputs/vcompare/v018/
( cd test-outputs/vcompare/v018 && \
  ../../../prolfquapp_docker.sh --image-version 0.1.8 prolfqua_yaml.sh -n robscale -y config.yaml -w cmp018 -s MAXQUANT -o . && \
  ../../../prolfquapp_docker.sh --image-version 0.1.8 prolfqua_dea.sh  -i . -d dataset.csv -y config.yaml -s MAXQUANT -o . -w cmp018 )

# 2.x: reuse the fixture's 2.x config and the prefixed software key
for V in 2.0.19 2.2.6; do
  d=test-outputs/vcompare/v${V//./}; mkdir -p "$d" && cp -r fixtures/maxquant_ionstar/* "$d/"
  ( cd "$d" && ../../../prolfquapp_docker.sh --image-version "$V" prolfqua_dea.sh \
      -i . -d dataset.csv -y config.yaml -s prolfquapp.MAXQUANT -o . -w "cmp${V//./}" )
done

Rscript scripts/compare_docker_versions.R \
  "0.1.8=$(find test-outputs/vcompare/v018  -name 'DE_*.xlsx' | head -1)" \
  "2.0.19=$(find test-outputs/vcompare/v2019 -name 'DE_*.xlsx' | head -1)" \
  "2.2.6=$(find test-outputs/vcompare/v226  -name 'DE_*.xlsx' | head -1)"
```

On ARM Mac these run under `linux/amd64` emulation (the older images are amd64-only).

### Result (maxquant_ionstar, contrast `e_vs_b`, 100 proteins)

| Pair | diff (log2FC) | p.value | FDR |
|------|---------------|---------|-----|
| 0.1.8 → 2.0.19 | 0.999916 (max\|Δ\| 0.011) | 0.9814 (0.23) | 0.9756 (0.08) |
| 2.0.19 → 2.2.6 | **1.000000 (0)** | **1.000000 (0)** | **1.000000 (0)** |
| 0.1.8 → 2.2.6 | 0.999916 (max\|Δ\| 0.011) | 0.9814 (0.23) | 0.9756 (0.08) |

**Interpretation:** `2.0.19 → 2.2.6` is **bit-identical** across all three metrics — the contrast-schema
refactor in `2.2.6` (prolfqua `1.6.3`: `modelName` = facade key, new `estimate_type` column, registry
consolidation) is purely structural and does not move fold-change or p-value estimates. The only real
change is `0.1.x → 2.x`: fold-changes stay essentially identical (r ≈ 0.9999) while p-value/FDR shift
modestly (r ≈ 0.98), reflecting genuine modelling/moderation evolution across the major-version jump.

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
| bgs_site_phospho | prolfquappPTMreaders | `inst/extdata/BGS_site/Report_test.tsv` | Real sites expanded deterministically to 2 control + 2 treated acquisitions |
| diann_wu345302 | local fixture payload | DIA-NN report, FASTA files, dataset, and config template | WU345302 facade matrix |

## Regenerating fixtures

```bash
make clean              # remove generated fixtures and outputs, keep diann_wu345302
make fixtures           # regenerate from source data
```

`make fixtures` regenerates the standard MaxQuant, MSstats, FP_TMT, FP_singlesite, and BGS_site fixtures from
`prolfquadata`, `prolfquappPTMreaders`, and the Zenodo PTM archive. The BGS source report contains one real acquisition;
the generator retains its sites and metadata and creates two deterministic control and two deterministic treated
acquisitions by scaling `PTM.Quantity`. This fixture tests reader and DEA plumbing, not biological reproducibility.
The target does not currently recreate the local `diann_wu345302` payload used by `make wu345302-facades`.

`make clean` removes the generated standard fixtures, logs, and development outputs, but intentionally keeps
`fixtures/diann_wu345302` because that local DIA-NN payload is not recreated by `make fixtures`.
