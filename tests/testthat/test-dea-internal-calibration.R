# test-dea-internal-calibration.R
#
# Verifies that setting `internal` in the YAML config triggers center_to_reference
# normalization in the DEA pipeline (CMD_DEA_V2.R).
#
# The reference protein is sp|O43707|ACTN4_HUMAN, reliably present in the
# fragpipe_ionstar fixture (MSstats format).
#
# After centering, center_to_reference subtracts the per-sample value of the
# reference protein from every protein — so the reference protein's own values
# in the transformedData assay must be exactly 0 (within floating-point tolerance).

REFERENCE_PROTEIN <- "sp|O43707|ACTN4_HUMAN"
REFERENCE_PROTEIN_2 <- "sp|P04406|G3P_HUMAN"

# ---------------------------------------------------------------------------
# Helper: copy fixture to a temp dir, patch config.yaml, run CMD_DEA_V2.R
# ---------------------------------------------------------------------------
run_dea_with_internal <- function(reference_proteins) {
  fixture_path <- file.path(FIXTURE_DIR, "fragpipe_ionstar")
  stopifnot(dir.exists(fixture_path))

  prefix <- if (length(reference_proteins) > 1) "dea_internal_multi_" else "dea_internal_"
  workdir <- file.path(
    normalizePath(file.path(FIXTURE_DIR, "..", "test-outputs"), mustWork = FALSE),
    paste0(prefix, format(Sys.time(), "%H%M%S"))
  )
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

  file.copy(
    list.files(fixture_path, full.names = TRUE, recursive = FALSE),
    workdir,
    recursive = TRUE
  )

  cfg_src <- file.path(workdir, "config.yaml")
  cfg_lines <- readLines(cfg_src)
  insert_after <- grep("^  transform:", cfg_lines)
  stopifnot(length(insert_after) == 1)
  internal_yaml <- paste0(
    "  internal: [",
    paste0("'", reference_proteins, "'", collapse = ", "),
    "]"
  )
  cfg_lines <- append(cfg_lines, internal_yaml, after = insert_after)
  writeLines(cfg_lines, cfg_src)

  args <- c(
    "--vanilla",
    get_cmd_dea_v2_path(),
    "-i", workdir,
    "-d", file.path(workdir, "dataset.csv"),
    "-y", cfg_src,
    "-s", "prolfquapp.MSSTATS",
    "-o", workdir,
    "-w", "test_internal"
  )

  cmd <- paste(
    "cd", shQuote(workdir), "&&",
    "Rscript", paste(shQuote(args), collapse = " ")
  )

  result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  exit_code <- attr(result, "status")
  if (!is.null(exit_code) && exit_code != 0L) {
    message("DEA output:\n", paste(result, collapse = "\n"))
  }

  cat("DEA outputs in:", workdir, "\n", file = stderr())
  list(
    workdir = workdir,
    output = result,
    exit_code = if (is.null(exit_code)) 0L else exit_code
  )
}


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

check_single_reference_zero <- function(se, reference_protein) {
  expect_true(
    reference_protein %in% rownames(se),
    label = paste("Reference protein present in SE rownames:", reference_protein)
  )
  assay_mat <- SummarizedExperiment::assay(se, "transformedData")
  ref_row <- assay_mat[reference_protein, ]
  non_na <- ref_row[!is.na(ref_row)]
  expect_gt(length(non_na), 0,
    label = paste("Reference protein has non-NA values in transformedData:", reference_protein))
  expect_true(
    all(abs(non_na) < 1e-6),
    label = paste(
      "Reference protein transformedData values are all ~0 after centering.",
      "Values:", paste(round(non_na, 8), collapse = ", ")
    )
  )
}

check_reference_median_zero <- function(se, reference_proteins) {
  # center_to_reference subtracts the per-sample MEDIAN of all reference proteins.
  # So the per-sample median of reference rows must be ~0; individual rows need not be.
  for (ref in reference_proteins) {
    expect_true(ref %in% rownames(se),
      label = paste("Reference protein present in SE rownames:", ref))
  }
  assay_mat <- SummarizedExperiment::assay(se, "transformedData")
  ref_mat <- assay_mat[reference_proteins, , drop = FALSE]
  per_sample_median <- apply(ref_mat, 2, median, na.rm = TRUE)
  non_na <- per_sample_median[!is.na(per_sample_median)]
  expect_gt(length(non_na), 0,
    label = "Reference proteins have non-NA values in transformedData")
  expect_true(
    all(abs(non_na) < 1e-6),
    label = paste(
      "Per-sample median of reference proteins is ~0 after centering.",
      "Values:", paste(round(non_na, 8), collapse = ", ")
    )
  )
}

test_that("internal calibration via CMD_DEA_V2.R: single reference protein centered to zero", {
  skip_on_cran()

  res <- run_dea_with_internal(REFERENCE_PROTEIN)

  expect_equal(res$exit_code, 0L,
    info = paste("DEA failed. Output:\n", paste(res$output, collapse = "\n")))

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$se_rds), 0, label = "SummarizedExperiment.rds produced")

  se <- readRDS(outputs$se_rds[1])
  expect_s4_class(se, "SummarizedExperiment")

  check_single_reference_zero(se, REFERENCE_PROTEIN)
})

test_that("internal calibration via CMD_DEA_V2.R: list of reference proteins - per-sample median centered to zero", {
  skip_on_cran()

  ref_proteins <- c(REFERENCE_PROTEIN, REFERENCE_PROTEIN_2)
  res <- run_dea_with_internal(ref_proteins)

  expect_equal(res$exit_code, 0L,
    info = paste("DEA failed. Output:\n", paste(res$output, collapse = "\n")))

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$se_rds), 0, label = "SummarizedExperiment.rds produced")

  se <- readRDS(outputs$se_rds[1])
  expect_s4_class(se, "SummarizedExperiment")

  check_reference_median_zero(se, ref_proteins)
})
