# CLAUDE.md — integration_test

## Purpose

Cross-package integration tests for the prolfquapp CLI pipelines (`prolfqua_dea.sh` / `CMD_DEA_V2.R` and `prolfqua_qc.sh` / `CMD_QUANT_QC.R`). Tests run the full pipeline end-to-end using small fixture datasets (~50–100 proteins) subsetted from real data.

## Build & Test

```bash
make help                    # Show all targets
make install                 # Reinstall prolfqua + prolfquapp + prolfquappPTMreaders from local source
make fixtures                # Generate fixture data from real datasets (one-time, ~3 min)
make test                    # Run all integration tests (~2-3 min)
make test-dea-maxquant       # Run a single test file
make save-references-docker  # Generate regression baselines from released Docker image
make test-dea-regression     # Compare dev output against Docker-generated baselines
```

## Regression References

The regression test (`test-dea-regression.R`) compares SummarizedExperiment outputs against saved baselines. **References must come from the released Docker image** (`prolfqua/prolfquapp:latest` on Docker Hub), not the local dev version.

```bash
make save-references-docker   # Run DEA in Docker, save reference SE.rds files
make test-dea-regression      # Compare local dev output against those references
```

On ARM Mac, Docker runs with `--platform linux/amd64` (Rosetta emulation).

The `fp_singlesite_phospho` fixture is skipped in Docker reference generation because the released image does not include `prolfquappPTMreaders`. The regression test skips gracefully when the reference is missing.

## Key Notes

- All DEA tests use `CMD_DEA_V2.R` (the legacy `CMD_DEA.R` was removed)
- QC tests use `CMD_QUANT_QC.R`
- Software naming: `CMD_DEA_V2.R` uses prefixed keys (`prolfquapp.MAXQUANT`), `CMD_QUANT_QC.R` uses unprefixed keys (`MAXQUANT`)
- Docker image: `prolfqua/prolfquapp:latest` — R libs at `/opt/r-libs-site/`, scripts at `/opt/r-libs-site/prolfquapp/application/`
