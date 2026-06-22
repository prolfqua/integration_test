test_that("DIANN protein reader auto-switches for a nested facade on real data", {
  skip_on_cran()
  skip_if_not_installed("limpa")

  fixture_path <- file.path(FIXTURE_DIR, "diann_wu345302")
  skip_if_not(dir.exists(fixture_path), "diann_wu345302 fixture not installed")

  output_base <- normalizePath(file.path(FIXTURE_DIR, "..", "test-outputs"), mustWork = FALSE)
  workdir <- file.path(
    output_base,
    paste0("dea_diann_auto_peptide_", format(Sys.time(), "%H%M%S"))
  )
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

  file.copy(
    list.files(fixture_path, full.names = TRUE, recursive = FALSE),
    workdir,
    recursive = TRUE
  )

  config_template <- file.path(workdir, "configs", "config_template.yaml")
  config_file <- file.path(workdir, "config.yaml")
  config <- readLines(config_template)
  config[grepl("^[[:space:]]*model:", config)] <- "  model: limpa_nested"
  writeLines(config, config_file)

  args <- c(
    "--vanilla",
    get_cmd_dea_v2_path(),
    "-i", file.path(workdir, "out-DIANN"),
    "-d", file.path(workdir, "dataset_saint.csv"),
    "-y", config_file,
    "-s", "prolfquapp.DIANN",
    "-o", workdir,
    "-w", "auto_peptide"
  )

  cmd <- paste(
    "cd", shQuote(workdir), "&&",
    "Rscript", paste(shQuote(args), collapse = " ")
  )
  message("Running DIANN auto-peptide DEA: ", cmd)

  result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  exit_code <- attr(result, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    message("DEA output:\n", paste(result, collapse = "\n"))
  }
  expect_equal(if (is.null(exit_code)) 0L else exit_code, 0L)

  outputs <- find_dea_outputs(workdir)
  expect_gt(length(outputs$se_rds), 0)
  expect_gt(length(outputs$yaml), 0)

  run_output <- paste(result, collapse = "\n")
  expect_match(run_output, "switching software from 'prolfquapp.DIANN'")
  expect_match(run_output, "prolfquapp.DIANN_PEPTIDE")

  minimal <- yaml::read_yaml(outputs$yaml[1])
  expect_equal(minimal$software, "prolfquapp.DIANN_PEPTIDE")

  log_files <- grep("prolfqua_.*\\.log$", outputs$all, value = TRUE)
  expect_gt(length(log_files), 0)
  log_text <- paste(unlist(lapply(log_files, readLines, warn = FALSE)), collapse = "\n")
  expect_match(log_text, "Processing done: prolfquapp.DIANN_PEPTIDE")
})
