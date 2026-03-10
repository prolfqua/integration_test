export NOT_CRAN=true

.DEFAULT_GOAL := help

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-25s\033[0m %s\n", $$1, $$2}'

LOGDIR := logs

$(LOGDIR):
	mkdir -p $(LOGDIR)

test: fixtures $(LOGDIR)  ## Run all integration tests
	Rscript -e "testthat::test_dir('tests/testthat')" 2>&1 | tee $(LOGDIR)/test-all.log

test-dea-maxquant: fixtures $(LOGDIR)  ## Test DEA with MaxQuant preprocessor
	Rscript -e "testthat::test_file('tests/testthat/test-dea-maxquant.R')" 2>&1 | tee $(LOGDIR)/test-dea-maxquant.log

test-dea-msstats: fixtures $(LOGDIR)  ## Test DEA with MSstats preprocessor
	Rscript -e "testthat::test_file('tests/testthat/test-dea-msstats.R')" 2>&1 | tee $(LOGDIR)/test-dea-msstats.log

test-dea-fp-tmt: fixtures $(LOGDIR)  ## Test DEA with FP_TMT preprocessor
	Rscript -e "testthat::test_file('tests/testthat/test-dea-fp-tmt.R')" 2>&1 | tee $(LOGDIR)/test-dea-fp-tmt.log

test-dea-fp-singlesite: fixtures $(LOGDIR)  ## Test DEA with FP_singlesite (phospho)
	Rscript -e "testthat::test_file('tests/testthat/test-dea-fp-singlesite.R')" 2>&1 | tee $(LOGDIR)/test-dea-fp-singlesite.log

test-qc-maxquant: fixtures $(LOGDIR)  ## Test QC pipeline with MaxQuant
	Rscript -e "testthat::test_file('tests/testthat/test-qc-maxquant.R')" 2>&1 | tee $(LOGDIR)/test-qc-maxquant.log

test-dea-regression: fixtures $(LOGDIR)  ## Test DEA outputs against saved references
	Rscript -e "testthat::test_file('tests/testthat/test-dea-regression.R')" 2>&1 | tee $(LOGDIR)/test-dea-regression.log

test-dea-internal: fixtures $(LOGDIR)  ## Test internal calibration (center_to_reference)
	Rscript -e "testthat::test_file('tests/testthat/test-dea-internal-calibration.R')" 2>&1 | tee $(LOGDIR)/test-dea-internal.log
	@echo "--- DEA output directories ---" && grep "DEA outputs in:" $(LOGDIR)/test-dea-internal.log || true

fixtures: fixtures/.stamp  ## Generate fixture data from real datasets (~3 min)

fixtures/.stamp:
	Rscript R/create_test_fixtures.R
	touch $@

save-references: fixtures  ## Save baseline SummarizedExperiment objects for regression
	Rscript -e "library(testthat); source('tests/testthat/helper-common.R'); source('tests/testthat/test-dea-regression.R'); save_references()"

install:  ## Reinstall prolfqua + prolfquapp + prolfquappPTMreaders from local source
	Rscript -e "devtools::install('../prolfqua', quick = TRUE, upgrade = 'never')"
	Rscript -e "devtools::install('../prolfquapp', quick = TRUE, upgrade = 'never')"
	Rscript -e "devtools::install('../prolfquappPTMreaders', quick = TRUE, upgrade = 'never')"

clean:  ## Remove all fixtures and start fresh
	rm -rf fixtures
	rm -rf tests/testthat/_snaps
	rm -rf $(LOGDIR)
	rm -rf test-outputs

clean-references:  ## Remove only regression references (keeps fixtures)
	rm -rf fixtures/reference

.PHONY: help test test-dea-maxquant test-dea-msstats test-dea-fp-tmt test-dea-fp-singlesite test-qc-maxquant test-dea-regression test-dea-internal save-references install clean clean-references
