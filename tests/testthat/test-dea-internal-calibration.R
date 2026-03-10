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

# ---------------------------------------------------------------------------
# Helper: copy fixture to a temp dir, patch config.yaml, run CMD_DEA_V2.R
# ---------------------------------------------------------------------------
run_dea_with_internal <- function(reference_protein) {
  fixture_path <- file.path(FIXTURE_DIR, "fragpipe_ionstar")
  stopifnot(dir.exists(fixture_path))

  workdir <- file.path(
    normalizePath(file.path(FIXTURE_DIR, "..", "test-outputs"), mustWork = FALSE),
    paste0("dea_internal_", format(Sys.time(), "%H%M%S"))
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
  cfg_lines <- append(cfg_lines,
    paste0("  internal: ['", reference_protein, "']"),
    after = insert_after)
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

test_that("internal calibration via CMD_DEA_V2.R: pipeline completes and centers reference protein to zero", {
  skip_on_cran()

  res <- run_dea_with_internal(REFERENCE_PROTEIN)

  expect_equal(res$exit_code, 0L,
    info = paste("DEA failed. Output:\n", paste(res$output, collapse = "\n")))

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$se_rds), 0, label = "SummarizedExperiment.rds produced")

  se <- readRDS(outputs$se_rds[1])
  expect_s4_class(se, "SummarizedExperiment")

  expect_true(
    REFERENCE_PROTEIN %in% rownames(se),
    label = paste("Reference protein present in SE rownames:", REFERENCE_PROTEIN)
  )

  # transformedData holds center_to_reference-adjusted values.
  # At protein level, each sample has one value for the reference protein.
  # Centering subtracts that value from itself → must be 0.
  assay_mat <- SummarizedExperiment::assay(se, "transformedData")
  ref_row <- assay_mat[REFERENCE_PROTEIN, ]
  non_na <- ref_row[!is.na(ref_row)]

  expect_gt(length(non_na), 0,
    label = "Reference protein has non-NA values in transformedData")

  expect_true(
    all(abs(non_na) < 1e-6),
    label = paste(
      "Reference protein transformedData values are all ~0 after centering.",
      "Values:", paste(round(non_na, 8), collapse = ", ")
    )
  )

})
