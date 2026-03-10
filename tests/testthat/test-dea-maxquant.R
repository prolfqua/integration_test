test_that("DEA pipeline runs with MAXQUANT software on IonStar fixture", {
  skip_on_cran()

  res <- run_dea(
    fixture_name = "maxquant_ionstar",
    software = "prolfquapp.MAXQUANT",
    workunit = "test_mq"
  )

  expect_equal(res$exit_code, 0L,
    info = paste("DEA failed. Output:\n", paste(res$output, collapse = "\n")))

  outputs <- find_dea_outputs(res$workdir)

  # Check that key output files exist
  expect_gt(length(outputs$html), 0, label = "HTML report(s) produced")
  expect_gt(length(outputs$xlsx), 0, label = "XLSX file(s) produced")
  expect_gt(length(outputs$se_rds), 0, label = "SummarizedExperiment.rds produced")
  expect_gt(length(outputs$rnk), 0, label = "RNK file(s) produced")
  expect_gt(length(outputs$parquet), 0, label = "Parquet file(s) produced")

  # Load SummarizedExperiment and check structure
  se <- readRDS(outputs$se_rds[1])
  expect_s4_class(se, "SummarizedExperiment")

  # Should have proteins (rows) and samples (columns)
  expect_gt(nrow(se), 0, label = "SE has proteins")
  expect_gt(ncol(se), 0, label = "SE has samples")

  # Contrast results stored in nested DataFrame columns (constrast_*)
  rd <- SummarizedExperiment::rowData(se)
  contrast_cols <- grep("^constrast_", colnames(rd), value = TRUE)
  expect_gt(length(contrast_cols), 0, label = "Contrast columns in rowData")

  # Check that the nested contrast DataFrame has fold-change (diff) and FDR
  first_contrast <- rd[[contrast_cols[1]]]
  expect_true("diff" %in% colnames(first_contrast), label = "diff column in contrast")
  expect_true("FDR" %in% colnames(first_contrast), label = "FDR column in contrast")

  # Should have non-NA fold-change values
  expect_gt(sum(!is.na(first_contrast$diff)), 0, label = "Non-NA diff values")

  # Clean up
  unlink(res$workdir, recursive = TRUE)
})
