# Integration Test Infrastructure for the prolfqua Ecosystem

## Context

You're planning a refactoring of `prolfqua` and need integration tests that catch breakage across the downstream package chain. Currently there are no cross-package integration tests — each package tests in isolation. The ecosystem has 5 interdependent R packages plus a Python/Snakemake pipeline, each in separate git repos under `~/projects/`.

**The problem:** A change to `prolfqua::LFQData` or any R6 class can silently break `prolfquapp`, `prophosqua`, and `prolfquappPTMreaders`. 
COMMENT: you did not include prolfquadata and prolfquabenchmark!!

There is no automated way to detect this before pushing.

## Decision: Standalone Integration Test Project (NOT a Monorepo)

**Create `~/projects/prolfqua-integration-tests/`** — a new git repo that orchestrates builds, checks, and integration tests across the package chain via Snakemake.

Comment: this is find with me. I still asked on the disc organziation, should we have a subfolder project/prolfqua/prolfqua, prolfquabenchmark, etc etc.prolfqua-integration-tests. Still this should be independent repos.



Why NOT a `prolfqua_projects/` monorepo with submodules:
- All packages are actively developed in their own repos and pushed independently
- Git submodules add daily friction (detached HEAD, `submodule update` required)
- The integration project needs to be a lightweight orchestrator, not own the code
- `config.yaml` pointing to sibling dirs (`../prolfqua`, etc.) is simpler and overridable for CI

## Package Dependency Graph

Comment: the dependency structure is wrong. prolfquabenchmark depends I think only on prolfqua.


```
Layer 1:  prolfqua                    (github.com/fgcz/prolfqua, branch: Modelling2R6)
             │
     ┌───────┴───────┐
Layer 2:  prolfquapp       prophosqua       (can install/check in parallel)
             │                  │
     └───────┬───────┘
Layer 3:  prolfquappPTMreaders              (depends on ALL above)


Also:     prolfquabenchmark                 (depends on prolfqua + prolfquapp + prolfquadata@gitlab.bfabric.org)
          ptm-pipeline (Python)             (runtime R dependency on prolfquapp + prophosqua)
```

## Test

run individual package test upon changes,
run all the sh files implemented in proflquapp against test data. Since the outputs are xlsx files, summarized experiments etc we need a method to validate that for our default setting they stay unchanged.



## Snakemake Workflow — Execution DAG

comment : Can you consult prolfquapPackagesBuildScript folder. maybe we should run all this against a pristine R lib? 

```
                    install_prolfqua
                   /        |        \
                  /         |         \
  install_prolfquapp  install_prophosqua  check_prolfqua prolfquappPTMreaders etc.
```

With `-j4`: Layer 2 installs + Layer 1 check run in parallel. Integration tests run in parallel after their deps.

### Key Snakemake Rules

Each package gets an `install_<pkg>` rule (runs `devtools::install()`) and a `check_<pkg>` rule (runs `devtools::check()` via `scripts/run_check.R`). Install rules depend on upstream installs. Check rules depend on their own install. Integration test rules depend on the relevant installs.

Convenience targets: `rule all`, `rule checks`, `rule integration_tests`, `rule installs`, `rule clean`.

### Makefile Wrapper

```makefile
JOBS ?= 4
all:           snakemake -j$(JOBS) all
checks:        snakemake -j$(JOBS) checks
integration:   snakemake -j$(JOBS) integration_tests
installs:      snakemake -j$(JOBS) installs
dry-run:       snakemake -n all
clean:         snakemake clean
```

## Integration Test Scripts

Comment: do not think so the most critical test is to ensure that all the scripts in prolfquapp do run and that the ptm piplines keeps running on reganerated intput data.


### 1. `smoke_test_prolfquapp.R` — Most Critical

comment: redundant to prophosqua.

Exercises the exact code path that breaks when prolfqua R6 interfaces change:
- `prolfqua::sim_lfq_data_peptide_config()` → `LFQData$new()` → `prolfquapp::ProteinAnnotation$new()` → `prolfquapp::DEAnalyse$new()` → `remove_cont_decoy()` → `aggregate()` → `transform_data()` → `build_model_linear_protein()` → `get_contrasts_linear_protein()` → `get_contrasts_merged_protein()`
- Verifies output structure: data.frame with columns `contrast`, `diff`, `statistic`, `p.value`, `FDR`

### 2. `verify_exports.R` — API Surface Regression

comment: unclear which package you are talking about. Interfaces are checked by running the dependent packages.

Checks that all R6 classes, public methods, and active bindings that downstream packages depend on still exist:
- 15+ critical prolfqua R6 classes (LFQData, Contrasts*, etc.)
- 20+ LFQData public methods/active bindings
- 8+ critical prolfquapp symbols (DEAnalyse, ProteinAnnotation, etc.)
- Fails immediately on any missing export with a clear message

comment: from now on I honestly not sure if you did not lost it completely. We are not in "Wasserfall model" you plan you do it. we proceed step by step and then we think if we need more.
So I am rather sure you can delte the rest for the moment, and focus on the above.


### 3. `smoke_test_fp.R` — Plugin Architecture

Tests the `prolfqua_preprocess_functions` plugin resolution, file discovery, and preprocessing with FragPipe sample data from `test_data/fragpipe/`.

### 4. `smoke_test_prophosqua.R` — PTM Analysis

Loads bundled prophosqua example data, exercises `prepare_site_protein_features()` and DPA/DPU utilities.

### 5. `smoke_test_ptmreaders.R` — Full 4-Package Chain

Loads all 4 packages, tests each PTM reader format (FP_multisite, FP_combined_STY, BGS_site), verifies function signature compatibility with `prolfquapp::preprocess_dummy()`.

### 6. `smoke_test_cli_dea.sh` — End-to-End CLI

Runs `CMD_DEA.R` with a minimal annotation and config, verifies output files are created.

## Test Data Strategy

| Source | Size | Action |
|--------|------|--------|
| prolfqua/data/*.rda | 2 MB | Keep in package (unit tests use it) |
| prolfquapp/inst/samples/FragPipe/ | 1 MB | Copy to test_data/fragpipe/ |
| prolfquapp/inst/samples/maxquant_txt/ | ~2 MB subset | Copy subset to test_data/maxquant/ |
| prolfquappPTMreaders/inst/extdata/ | 200 KB | Copy to test_data/ptm_*/ |
| prophosqua/data/*.rda | 11 MB | Keep in package (loaded via `data()`) |
| prolfquapp/inst/application/DIANN/ | 1.8 GB | **Never in CI** — stay local |
| ptm-pipeline/test_data/ | 12 GB | **Never in CI** — stay local |
| New: diann_tiny/ | ~1 MB | Source from ProteoBench test data (small DIANN output) |

**Budget: < 50 MB total in test_data/**, all git-tracked.

## GitHub Actions CI

File: `.github/workflows/integration-tests.yml`

- **Triggers:** push/PR to main, nightly schedule (weekdays 03:00 UTC), manual dispatch with branch override inputs, `repository_dispatch` from upstream repos
- **Steps:** checkout → setup R 4.5.2 → setup Python 3.12 → install Snakemake → cache R library → install system deps → install R deps (BiocManager for vsn, limma, etc.) → clone ecosystem packages (shallow, configurable branches) → `snakemake -j2 all` → upload logs on failure
- **Caching:** R library cached keyed on dependency list hash

Optional: Add `trigger-integration.yml` to each upstream repo to fire `repository_dispatch` on push to main branch.

## Implementation Phases

### Phase 1: Bootstrap
1. `mkdir ~/projects/prolfqua-integration-tests && cd $_ && git init`
2. Create `config.yaml`, `Snakefile`, `Makefile`, `CLAUDE.md`, `.gitignore`
3. Create `scripts/run_check.R`

### Phase 2: Test Data
1. Create `test_data/` subdirs, copy small files from existing packages
2. Create annotation CSVs for each test scenario
3. Write `test_data/README.md` with provenance

### Phase 3: Integration Scripts
1. `scripts/smoke_test_prolfquapp.R` (highest priority — DEAnalyse pipeline)
2. `scripts/verify_exports.R` (API regression)
3. `scripts/smoke_test_fp.R`
4. `scripts/smoke_test_prophosqua.R`
5. `scripts/smoke_test_ptmreaders.R`
6. `scripts/smoke_test_cli_dea.sh`

### Phase 4: Local Validation
1. `make dry-run` → verify DAG
2. `make all` → fix issues
3. Iterate until green

### Phase 5: CI
1. Create `.github/workflows/integration-tests.yml`
2. Push to GitHub, verify
3. (Later) Add upstream dispatch triggers

## Critical Files to Reference During Implementation

- `prolfquapp/R/R6_DEAnalyse.R` — DEAnalyse class, primary consumer of prolfqua API
- `prolfquapp/R/preprocess_software.R` — Plugin architecture, `prolfqua_preprocess_functions` list
- `prolfqua/R/LFQData.R` — LFQData R6 class, core API surface to protect
- `prolfquappPTMreaders/R/*.R` — Reader functions with signature contracts
- `prolfqua/.github/workflows/r.yaml` — Existing CI pattern to follow
- `proLFQuaPackageBuildScripts/buildProlfqua.sh` — Existing sequential build script (being replaced)

## Verification

After implementation, verify end-to-end:
1. `cd ~/projects/prolfqua-integration-tests`
2. `make dry-run` — should show the full DAG with no errors
3. `make installs` — all 4 packages install successfully
4. `make checks` — all 4 R CMD checks pass (warnings allowed initially)
5. `make integration` — all smoke tests pass
6. `make all` — full suite green
7. Intentionally break something in prolfqua (e.g., rename an LFQData method) and confirm the integration tests catch it
