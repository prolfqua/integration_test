# History

## 2026-02-21 — Build infrastructure

- Synchronized Makefiles across all R packages (common targets: test, check, document, lint, format, renv-init/restore/snapshot)
- Added `Suggests` to DESCRIPTION in packages that were missing it (prolfquappPTMreaders, prolfquasaint, prolfquabenchmark)
- Added dev tools (devtools, roxygen2, covr, lintr, pkgdown, rcmdcheck) to Suggests in all 6 R packages — required because renv isolates from system libraries
- Documented hybrid renv strategy in CLAUDE.md: top-level renv for cross-package dev, per-package renv for isolated checks
- Fixed CLAUDE.md rule: NAMESPACE is generated (don't edit), DESCRIPTION is edited directly
- Initialized renv in prolfqua, confirmed `make check` runs against it
