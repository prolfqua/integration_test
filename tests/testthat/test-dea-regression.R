# test-dea-regression.R
#
# Regression tests: compare SummarizedExperiment outputs against saved references.
# For every nested DataFrame in rowData and for the assay matrix, we correlate
# all numeric columns against the reference and expect r > 0.999.
#
# References are created by: make save-references

# FIXTURE_DIR is an absolute path (from helper-common.R), so this is stable
REFERENCE_DIR <- normalizePath(file.path(FIXTURE_DIR, "..", "reference"),
                               mustWork = FALSE)

# Helper to save reference SummarizedExperiment objects
save_references <- function() {
  dir.create(REFERENCE_DIR, recursive = TRUE, showWarnings = FALSE)

  fixtures <- list(
    maxquant = list(name = "maxquant_ionstar", software = "prolfquapp.MAXQUANT",
                    dataset = "dataset.csv", workunit = "ref_mq"),
    msstats = list(name = "fragpipe_ionstar", software = "prolfquapp.MSSTATS",
                   dataset = "dataset.csv", workunit = "ref_msstats"),
    fp_tmt = list(name = "fp_tmt_total", software = "prolfquapp.FP_TMT",
                  dataset = "dataset_with_contrasts.tsv", workunit = "ref_fp_tmt"),
    fp_singlesite = list(name = "fp_singlesite_phospho",
                         software = "prolfquappPTMreaders.FP_singlesite",
                         dataset = "dataset_with_contrasts.tsv",
                         workunit = "ref_fp_singlesite")
  )

  for (key in names(fixtures)) {
    f <- fixtures[[key]]
    message("Running DEA for reference: ", key)
    res <- run_dea(f$name, f$software, f$dataset, workunit = f$workunit)
    if (res$exit_code != 0L) {
      warning("DEA failed for ", key, ". Skipping reference save.")
      next
    }
    outputs <- find_dea_outputs(res$workdir)
    if (length(outputs$se_rds) > 0) {
      out_file <- file.path(REFERENCE_DIR, paste0(key, "_SE.rds"))
      file.copy(outputs$se_rds[1], out_file, overwrite = TRUE)
      message("  Saved: ", out_file)
    }
    unlink(res$workdir, recursive = TRUE)
  }
  message("All references saved to: ", REFERENCE_DIR)
}


# ---------------------------------------------------------------------------
# Extract all numeric vectors from a nested rowData DataFrame
# Returns a named list of numeric vectors (name = "parent.column")
# ---------------------------------------------------------------------------
extract_numeric_vectors <- function(rd) {
  result <- list()
  for (col in colnames(rd)) {
    val <- rd[[col]]
    if (is.data.frame(val) || inherits(val, "DataFrame")) {
      # Nested DataFrame (e.g. constrast_*, stats_*_wide)
      for (subcol in colnames(val)) {
        v <- tryCatch(as.numeric(val[[subcol]]), warning = function(w) NULL)
        if (!is.null(v) && sum(!is.na(v)) > 0) {
          result[[paste0(col, ".", subcol)]] <- v
        }
      }
    } else {
      v <- tryCatch(as.numeric(val), warning = function(w) NULL)
      if (!is.null(v) && sum(!is.na(v)) > 0) {
        result[[col]] <- v
      }
    }
  }
  result
}


# ---------------------------------------------------------------------------
# Compare a new SE against a reference SE
# Checks: dimensions, column names, correlation of all numeric data >= min_cor
# ---------------------------------------------------------------------------
compare_se <- function(new_se, ref_se, min_cor = 0.999) {
  # Same dimensions
  expect_equal(nrow(new_se), nrow(ref_se), label = "Number of rows")
  expect_equal(ncol(new_se), ncol(ref_se), label = "Number of columns")

  # --- Assay matrix correlation ---
  new_mat <- as.numeric(SummarizedExperiment::assay(new_se))
  ref_mat <- as.numeric(SummarizedExperiment::assay(ref_se))
  both_valid <- !is.na(new_mat) & !is.na(ref_mat)
  if (sum(both_valid) > 2) {
    r <- cor(new_mat[both_valid], ref_mat[both_valid])
    expect_gte(r, min_cor,
      label = sprintf("Assay matrix correlation (r=%.6f)", r))
  }

  # --- rowData: correlate every numeric column ---
  new_rd <- SummarizedExperiment::rowData(new_se)
  ref_rd <- SummarizedExperiment::rowData(ref_se)

  expect_equal(sort(colnames(new_rd)), sort(colnames(ref_rd)),
    label = "rowData column names")

  new_vecs <- extract_numeric_vectors(new_rd)
  ref_vecs <- extract_numeric_vectors(ref_rd)

  # Check all reference numeric columns exist in new
  shared <- intersect(names(new_vecs), names(ref_vecs))
  expect_gt(length(shared), 0, label = "Shared numeric columns")

  for (col in shared) {
    nv <- new_vecs[[col]]
    rv <- ref_vecs[[col]]
    both_valid <- !is.na(nv) & !is.na(rv)
    if (sum(both_valid) > 2) {
      # Skip constant columns (sd = 0 → cor undefined)
      if (sd(rv[both_valid]) == 0 && sd(nv[both_valid]) == 0) next
      r <- cor(nv[both_valid], rv[both_valid])
      expect_gte(r, min_cor,
        label = sprintf("%s correlation (r=%.6f)", col, r))
    }
  }
}


test_that("MaxQuant DEA output matches reference", {
  skip_on_cran()
  ref_file <- file.path(REFERENCE_DIR, "maxquant_SE.rds")
  skip_if_not(file.exists(ref_file), "Reference not found. Run: make save-references")

  res <- run_dea("maxquant_ionstar", "prolfquapp.MAXQUANT", workunit = "reg_mq")
  expect_equal(res$exit_code, 0L)

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$se_rds), 0)

  new_se <- readRDS(outputs$se_rds[1])
  ref_se <- readRDS(ref_file)
  compare_se(new_se, ref_se)

  unlink(res$workdir, recursive = TRUE)
})


test_that("MSstats DEA output matches reference", {
  skip_on_cran()
  ref_file <- file.path(REFERENCE_DIR, "msstats_SE.rds")
  skip_if_not(file.exists(ref_file), "Reference not found. Run: make save-references")

  res <- run_dea("fragpipe_ionstar", "prolfquapp.MSSTATS", workunit = "reg_msstats")
  expect_equal(res$exit_code, 0L)

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$se_rds), 0)

  new_se <- readRDS(outputs$se_rds[1])
  ref_se <- readRDS(ref_file)
  compare_se(new_se, ref_se)

  unlink(res$workdir, recursive = TRUE)
})


test_that("FP_TMT DEA output matches reference", {
  skip_on_cran()
  ref_file <- file.path(REFERENCE_DIR, "fp_tmt_SE.rds")
  skip_if_not(file.exists(ref_file), "Reference not found. Run: make save-references")

  res <- run_dea("fp_tmt_total", "prolfquapp.FP_TMT",
                 dataset_file = "dataset_with_contrasts.tsv",
                 workunit = "reg_fp_tmt")
  expect_equal(res$exit_code, 0L)

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$se_rds), 0)

  new_se <- readRDS(outputs$se_rds[1])
  ref_se <- readRDS(ref_file)
  compare_se(new_se, ref_se)

  unlink(res$workdir, recursive = TRUE)
})


test_that("FP_singlesite DEA output matches reference", {
  skip_on_cran()
  skip_if_not_installed("prolfquappPTMreaders")
  ref_file <- file.path(REFERENCE_DIR, "fp_singlesite_SE.rds")
  skip_if_not(file.exists(ref_file), "Reference not found. Run: make save-references")

  res <- run_dea("fp_singlesite_phospho", "prolfquappPTMreaders.FP_singlesite",
                 dataset_file = "dataset_with_contrasts.tsv",
                 workunit = "reg_fp_singlesite")
  expect_equal(res$exit_code, 0L)

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$se_rds), 0)

  new_se <- readRDS(outputs$se_rds[1])
  ref_se <- readRDS(ref_file)
  compare_se(new_se, ref_se)

  unlink(res$workdir, recursive = TRUE)
})
