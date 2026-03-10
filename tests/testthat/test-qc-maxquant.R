test_that("QC pipeline runs with MAXQUANT software on IonStar fixture", {
  skip_on_cran()

  # Note: CMD_QUANT_QC.R uses prolfqua_preprocess_functions (unprefixed keys),
  # unlike CMD_DEA.R which uses get_procfuncs() (prefixed keys).
  # So QC needs plain "MAXQUANT", not "prolfquapp.MAXQUANT".
  res <- run_qc(
    fixture_name = "maxquant_ionstar",
    software = "MAXQUANT",
    workunit = "test_qc_mq"
  )

  expect_equal(res$exit_code, 0L,
    info = paste("QC failed. Output:\n", paste(res$output, collapse = "\n")))

  outputs <- find_qc_outputs(res$workdir)

  # QC produces HTML reports and XLSX
  expect_gt(length(outputs$html), 0, label = "QC HTML report(s) produced")
  expect_gt(length(outputs$xlsx), 0, label = "QC XLSX file(s) produced")

  unlink(res$workdir, recursive = TRUE)
})
