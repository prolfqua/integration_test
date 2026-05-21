export NOT_CRAN=true

.DEFAULT_GOAL := help

help:  ## Show this help
	@grep -E '^[a-zA-Z0-9_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-25s\033[0m %s\n", $$1, $$2}'

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

DOCKER_SH := $(CURDIR)/prolfquapp_docker.sh
REFDIR := reference

# Copy prolfquapp_docker.sh from the installed package (or Docker image) if not present
$(DOCKER_SH):
	R --vanilla -e "prolfquapp::copy_docker_script(workdir = '$(CURDIR)')"
	chmod +x $@

# Helper: run DEA in Docker for one fixture, extract SummarizedExperiment.rds
# prolfquapp_docker.sh mounts cwd at /work, so we set up a workdir with fixture data + output.
# Usage: $(call docker_dea,fixture_name,dataset_file,software,ref_key)
define docker_dea
	@echo "--- Generating reference: $(4) from fixture $(1) ---"
	rm -rf test-outputs/docker_$(4)
	mkdir -p test-outputs/docker_$(4)
	cp -r $(CURDIR)/fixtures/$(1)/* test-outputs/docker_$(4)/
	cd test-outputs/docker_$(4) && $(DOCKER_SH) prolfqua_dea.sh \
		-i . -d $(2) -y config.yaml \
		-s $(3) -o . -w ref_$(4)
	@se=$$(find test-outputs/docker_$(4) -name 'SummarizedExperiment.rds' | head -1); \
	if [ -z "$$se" ]; then echo "ERROR: No SummarizedExperiment.rds found for $(4)"; exit 1; fi; \
	cp "$$se" $(REFDIR)/$(4)_SE.rds; \
	echo "  Saved: $(REFDIR)/$(4)_SE.rds"
endef

save-references-docker: fixtures $(DOCKER_SH)  ## Generate regression references from released Docker image
	mkdir -p $(REFDIR)
	$(call docker_dea,maxquant_ionstar,dataset.csv,prolfquapp.MAXQUANT,maxquant)
	$(call docker_dea,fragpipe_ionstar,dataset.csv,prolfquapp.MSSTATS,msstats)
	$(call docker_dea,fp_tmt_total,dataset_with_contrasts.tsv,prolfquapp.FP_TMT,fp_tmt)
	@echo "=== References saved to $(REFDIR)/ (fp_singlesite skipped - not in Docker image) ==="

save-references: fixtures  ## Save baseline SE objects using local dev install (for re-baselining)
	Rscript -e "library(testthat); source('tests/testthat/helper-common.R'); source('tests/testthat/test-dea-regression.R'); save_references()"

install:  ## Reinstall prolfqua + prolfquapp + prolfquappPTMreaders from local source
	Rscript -e "devtools::install('../prolfqua', upgrade = 'never')"
	Rscript -e "devtools::install('../prolfquapp', build_vignettes = TRUE, upgrade = 'never')"
	Rscript -e "devtools::install('../prolfquappPTMreaders', upgrade = 'never')"

COMPARE_RMD := R/compare_regression.Rmd

# Fixture-to-reference mapping for compare-regression
COMPARE_FIXTURES := maxquant:maxquant_ionstar msstats:fragpipe_ionstar fp_tmt:fp_tmt_total

compare-regression:  ## Render visual comparison reports (ref vs dev) for each fixture
	@for entry in $(COMPARE_FIXTURES); do \
		key=$${entry%%:*}; fixture=$${entry##*:}; \
		ref=$(REFDIR)/$${key}_SE.rds; \
		if [ ! -f "$$ref" ]; then echo "SKIP $$key: no reference (run make save-references-docker)"; continue; fi; \
		test_se=$$(find test-outputs -path "*/dea_$${fixture}_*/**/SummarizedExperiment.rds" -newer "$$ref" 2>/dev/null | sort | tail -1); \
		if [ -z "$$test_se" ]; then \
			test_se=$$(find test-outputs -path "*/dea_$${fixture}_*/*/SummarizedExperiment.rds" 2>/dev/null | sort | tail -1); \
		fi; \
		if [ -z "$$test_se" ]; then echo "SKIP $$key: no test output (run make test first)"; continue; fi; \
		echo "--- Rendering comparison: $$key ---"; \
		echo "  ref:  $$ref"; \
		echo "  test: $$test_se"; \
		Rscript -e "rmarkdown::render('$(CURDIR)/$(COMPARE_RMD)', \
			params = list(ref_se_path = '$(CURDIR)/$$ref', \
			              test_se_path = '$(CURDIR)/$$test_se', \
			              title = 'Regression: $${key}'), \
			output_file = '$(CURDIR)/test-outputs/compare_$${key}.html')"; \
	done
	@echo "=== Comparison reports in test-outputs/compare_*.html ==="

wu345302-facades:  ## Run all registered WU345302 facade models and summarize correlations
	bash scripts/run_wu345302_facades.sh

clean:  ## Remove fixtures, logs, and dev test outputs (keeps Docker references)
	rm -rf fixtures
	rm -rf tests/testthat/_snaps
	rm -rf $(LOGDIR)
	find test-outputs -maxdepth 1 -mindepth 1 ! -name 'docker_*' -exec rm -rf {} + 2>/dev/null || true

clean-references:  ## Remove only regression references (keeps fixtures)
	rm -rf $(REFDIR)

.PHONY: help test test-dea-maxquant test-dea-msstats test-dea-fp-tmt test-dea-fp-singlesite test-qc-maxquant test-dea-regression test-dea-internal save-references save-references-docker compare-regression wu345302-facades install clean clean-references
