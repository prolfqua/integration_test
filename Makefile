export NOT_CRAN=true

.DEFAULT_GOAL := help

help:  ## Show this help
	@grep -E '^[a-zA-Z0-9_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-25s\033[0m %s\n", $$1, $$2}'

LOGDIR := logs

$(LOGDIR):
	mkdir -p $(LOGDIR)

# --- Functional / feature tests (run the locally installed prolfquapp) -------

test: fixtures $(LOGDIR)  ## Run all integration tests (local install)
	Rscript -e "testthat::test_dir('tests/testthat')" 2>&1 | tee $(LOGDIR)/test-all.log

test-dea-maxquant: fixtures $(LOGDIR)  ## Test DEA with MaxQuant preprocessor
	Rscript -e "testthat::test_file('tests/testthat/test-dea-maxquant.R')" 2>&1 | tee $(LOGDIR)/test-dea-maxquant.log

test-dea-msstats: fixtures $(LOGDIR)  ## Test DEA with MSstats preprocessor
	Rscript -e "testthat::test_file('tests/testthat/test-dea-msstats.R')" 2>&1 | tee $(LOGDIR)/test-dea-msstats.log

test-dea-fp-tmt: fixtures $(LOGDIR)  ## Test DEA with FP_TMT preprocessor
	Rscript -e "testthat::test_file('tests/testthat/test-dea-fp-tmt.R')" 2>&1 | tee $(LOGDIR)/test-dea-fp-tmt.log

test-dea-fp-singlesite: fixtures $(LOGDIR)  ## Test DEA with FP_singlesite (phospho)
	Rscript -e "testthat::test_file('tests/testthat/test-dea-fp-singlesite.R')" 2>&1 | tee $(LOGDIR)/test-dea-fp-singlesite.log

test-dea-bgs-site: fixtures $(LOGDIR)  ## Test DEA with Spectronaut BGS_site (phospho)
	Rscript -e "testthat::test_file('tests/testthat/test-dea-bgs-site.R')" 2>&1 | tee $(LOGDIR)/test-dea-bgs-site.log

test-qc-maxquant: fixtures $(LOGDIR)  ## Test QC pipeline with MaxQuant
	Rscript -e "testthat::test_file('tests/testthat/test-qc-maxquant.R')" 2>&1 | tee $(LOGDIR)/test-qc-maxquant.log

test-dea-diann-auto-peptide: fixtures $(LOGDIR)  ## Test nested facade auto-selects DIANN_PEPTIDE on real DIA-NN fixture
	Rscript -e "testthat::test_file('tests/testthat/test-dea-diann-auto-peptide.R')" 2>&1 | tee $(LOGDIR)/test-dea-diann-auto-peptide.log

test-dea-internal: fixtures $(LOGDIR)  ## Test internal calibration (center_to_reference)
	Rscript -e "testthat::test_file('tests/testthat/test-dea-internal-calibration.R')" 2>&1 | tee $(LOGDIR)/test-dea-internal.log
	@echo "--- DEA output directories ---" && grep "DEA outputs in:" $(LOGDIR)/test-dea-internal.log || true

# --- Fixtures ----------------------------------------------------------------

fixtures: fixtures/.stamp  ## Generate fixture data from real datasets (~3 min)

fixtures/.stamp: R/create_test_fixtures.R
	Rscript R/create_test_fixtures.R
	touch $@

install:  ## Reinstall prolfqua + prolfquapp + prolfquappPTMreaders from local source
	Rscript -e "devtools::install('../prolfqua', upgrade = 'never')"
	Rscript -e "devtools::install('../prolfquapp', build_vignettes = TRUE, upgrade = 'never')"
	Rscript -e "devtools::install('../prolfquappPTMreaders', upgrade = 'never')"

# --- Docker regression (released images only) --------------------------------
# Runs a fixture through released Docker tags / models and correlates the
# protein-level fold-change / p.value / FDR. See README "Docker regression".

regression: fixtures  ## Docker model-matrix regression (maxquant_ionstar)
	bash scripts/docker_regression.sh maxquant_ionstar prolfquapp.MAXQUANT dataset.csv 0.95

# --- WU345302 facade matrix --------------------------------------------------

wu345302-facades:  ## Run all registered WU345302 facade models and summarize correlations
	bash scripts/run_wu345302_facades.sh

# --- Cleanup -----------------------------------------------------------------

clean:  ## Remove generated fixtures, logs, and test outputs (keeps WU345302 fixture)
	rm -rf fixtures/.stamp
	rm -rf fixtures/maxquant_ionstar
	rm -rf fixtures/fragpipe_ionstar
	rm -rf fixtures/fp_tmt_total
	rm -rf fixtures/fp_singlesite_phospho
	rm -rf fixtures/bgs_site_phospho
	rm -rf tests/testthat/_snaps
	rm -rf $(LOGDIR)
	rm -rf test-outputs

.PHONY: help test fixtures install regression wu345302-facades clean
.PHONY: test-dea-maxquant test-dea-msstats test-dea-fp-tmt test-dea-fp-singlesite
.PHONY: test-dea-bgs-site test-qc-maxquant test-dea-diann-auto-peptide test-dea-internal
