test_that("DEA pipeline runs with FP_multisite on phospho fixture", {
  skip_on_cran()
  skip_if_not_installed("prolfquappPTMreaders")

  res <- run_dea(
    fixture_name = "fp_multisite_phospho",
    software = "prolfquappPTMreaders.FP_multisite",
    workunit = "test_fp_multisite"
  )

  expect_equal(
    res$exit_code,
    0L,
    info = paste("DEA failed. Output:\n", paste(res$output, collapse = "\n"))
  )

  outputs <- find_dea_outputs(res$workdir)
  expect_gt(length(outputs$html), 0, label = "HTML report(s) produced")
  expect_gt(length(outputs$xlsx), 0, label = "XLSX file(s) produced")
  expect_gt(
    length(outputs$se_rds),
    0,
    label = "SummarizedExperiment.rds produced"
  )

  se <- readRDS(outputs$se_rds[1])
  expect_s4_class(se, "SummarizedExperiment")
  expect_gt(nrow(se), 0, label = "SE has multisites")
  expect_equal(ncol(se), 8L, label = "SE has all eight multisite samples")

  row_data <- SummarizedExperiment::rowData(se)
  contrast_cols <- grep("^constrast_", colnames(row_data), value = TRUE)
  expect_gt(length(contrast_cols), 0, label = "Contrast columns in rowData")

  first_contrast <- row_data[[contrast_cols[1]]]
  expect_equal(
    setdiff(c("site", "diff", "FDR"), colnames(first_contrast)),
    character()
  )
  expect_complete_protein_annotation(first_contrast)
  expect_gt(sum(!is.na(first_contrast$diff)), 0, label = "Non-NA diff values")
})
