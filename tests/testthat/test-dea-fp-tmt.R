test_that("DEA pipeline runs with FP_TMT software on total proteome fixture", {
  skip_on_cran()

  res <- run_dea(
    fixture_name = "fp_tmt_total",
    software = "prolfquapp.FP_TMT",
    dataset_file = "dataset_with_contrasts.tsv",
    workunit = "test_fp_tmt"
  )

  expect_equal(res$exit_code, 0L,
    info = paste("DEA failed. Output:\n", paste(res$output, collapse = "\n")))

  outputs <- find_dea_outputs(res$workdir)

  expect_gt(length(outputs$html), 0, label = "HTML report(s) produced")
  expect_gt(length(outputs$xlsx), 0, label = "XLSX file(s) produced")
  expect_gt(length(outputs$se_rds), 0, label = "SummarizedExperiment.rds produced")
  expect_gt(length(outputs$rnk), 0, label = "RNK file(s) produced")

  se <- readRDS(outputs$se_rds[1])
  expect_s4_class(se, "SummarizedExperiment")
  expect_gt(nrow(se), 0, label = "SE has proteins")
  expect_gt(ncol(se), 0, label = "SE has samples")

  # This fixture has 4 complex contrasts (KO vs WT overall + per-timepoint)
  rd <- SummarizedExperiment::rowData(se)
  contrast_cols <- grep("^constrast_", colnames(rd), value = TRUE)
  expect_gte(length(contrast_cols), 4, label = "At least 4 contrast columns")

  first_contrast <- rd[[contrast_cols[1]]]
  expect_true("diff" %in% colnames(first_contrast))
  expect_true("FDR" %in% colnames(first_contrast))
  expect_gt(sum(!is.na(first_contrast$diff)), 0, label = "Non-NA diff values")

  # outputs persist in test-outputs/ for inspection
})
