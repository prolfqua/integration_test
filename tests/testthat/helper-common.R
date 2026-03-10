# helper-common.R
# Shared utilities for integration tests
# Loaded automatically by testthat before test files

library(prolfquapp)
library(prolfqua)

# Locate the fixture directory relative to this helper file.
# This works regardless of working directory (testthat, make, source()).
find_fixture_dir <- function() {
  # This file is at integration_test/tests/testthat/helper-common.R
  # Fixtures are at   integration_test/fixtures/
  helper_dir <- normalizePath(file.path("tests", "testthat"), mustWork = FALSE)
  if (!dir.exists(helper_dir)) {
    # testthat sets cwd to tests/testthat when running tests
    helper_dir <- normalizePath(".", mustWork = FALSE)
  }
  fixture_dir <- normalizePath(file.path(helper_dir, "..", "..", "fixtures"),
                               mustWork = FALSE)
  if (dir.exists(fixture_dir)) return(fixture_dir)
  stop("Cannot find fixtures directory at: ", fixture_dir,
       "\nRun 'make fixtures' from integration_test/")
}

FIXTURE_DIR <- find_fixture_dir()

# Get the path to CMD_DEA.R from the installed prolfquapp package
get_cmd_dea_path <- function() {
  pkg_path <- system.file(package = "prolfquapp")
  script <- file.path(pkg_path, "application", "CMD_DEA.R")
  stopifnot(file.exists(script))
  script
}

# Get the path to CMD_DEA_V2.R from the installed prolfquapp package
get_cmd_dea_v2_path <- function() {
  pkg_path <- system.file(package = "prolfquapp")
  script <- file.path(pkg_path, "application", "CMD_DEA_V2.R")
  stopifnot(file.exists(script))
  script
}

# Get the path to CMD_QUANT_QC.R from the installed prolfquapp package
get_cmd_qc_path <- function() {
  pkg_path <- system.file(package = "prolfquapp")
  script <- file.path(pkg_path, "application", "CMD_QUANT_QC.R")
  stopifnot(file.exists(script))
  script
}

# Run CMD_DEA.R in a temp directory with the given fixture
# Returns the output directory path
run_dea <- function(fixture_name,
                    software,
                    dataset_file = "dataset.csv",
                    config_file = "config.yaml",
                    workunit = NULL) {
  fixture_path <- file.path(FIXTURE_DIR, fixture_name)
  stopifnot(dir.exists(fixture_path))

  # Create a temporary working directory
  workdir <- file.path(tempdir(), paste0("dea_", fixture_name, "_", format(Sys.time(), "%H%M%S")))
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

  # Copy fixture contents to workdir
  file.copy(
    list.files(fixture_path, full.names = TRUE, recursive = FALSE),
    workdir,
    recursive = TRUE
  )

  # Build command args
  args <- c(
    "--vanilla",
    get_cmd_dea_path(),
    "-i", workdir,
    "-d", file.path(workdir, dataset_file),
    "-y", file.path(workdir, config_file),
    "-s", software,
    "-o", workdir
  )
  if (!is.null(workunit)) {
    args <- c(args, "-w", workunit)
  }

  # CMD_DEA.R's copy_DEA_Files() copies Rmd to cwd, so run from workdir
  cmd <- paste(
    "cd", shQuote(workdir), "&&",
    "Rscript", paste(shQuote(args), collapse = " ")
  )
  message("Running DEA: ", cmd)

  result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  exit_code <- attr(result, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    message("DEA output:\n", paste(result, collapse = "\n"))
  }

  list(
    workdir = workdir,
    output = result,
    exit_code = if (is.null(exit_code)) 0L else exit_code
  )
}

# Run CMD_QUANT_QC.R in a temp directory with the given fixture
run_qc <- function(fixture_name,
                   software,
                   dataset_file = "dataset.csv",
                   config_file = "config.yaml",
                   workunit = "test_qc") {
  fixture_path <- file.path(FIXTURE_DIR, fixture_name)
  stopifnot(dir.exists(fixture_path))

  workdir <- file.path(tempdir(), paste0("qc_", fixture_name, "_", format(Sys.time(), "%H%M%S")))
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

  file.copy(
    list.files(fixture_path, full.names = TRUE, recursive = FALSE),
    workdir,
    recursive = TRUE
  )

  args <- c(
    "--vanilla",
    get_cmd_qc_path(),
    "-i", workdir,
    "-d", file.path(workdir, dataset_file),
    "-y", file.path(workdir, config_file),
    "-s", software,
    "-o", file.path(workdir, "qc_output"),
    "-w", workunit
  )

  # CMD_QUANT_QC.R may copy files to cwd, so run from workdir
  cmd <- paste(
    "cd", shQuote(workdir), "&&",
    "Rscript", paste(shQuote(args), collapse = " ")
  )
  message("Running QC: ", cmd)

  result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  exit_code <- attr(result, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    message("QC output:\n", paste(result, collapse = "\n"))
  }

  list(
    workdir = workdir,
    output = result,
    exit_code = if (is.null(exit_code)) 0L else exit_code
  )
}

# Find output files produced by DEA in the workdir
find_dea_outputs <- function(workdir) {
  all_files <- list.files(workdir, recursive = TRUE, full.names = TRUE)
  list(
    html = grep("\\.html$", all_files, value = TRUE),
    xlsx = grep("\\.xlsx$", all_files, value = TRUE),
    rnk = grep("\\.rnk$", all_files, value = TRUE),
    parquet = grep("\\.parquet$", all_files, value = TRUE),
    se_rds = grep("SummarizedExperiment\\.rds$", all_files, value = TRUE),
    yaml = grep("minimal\\.yaml$", all_files, value = TRUE),
    all = all_files
  )
}

# Find output files produced by QC in the workdir
find_qc_outputs <- function(workdir) {
  all_files <- list.files(workdir, recursive = TRUE, full.names = TRUE)
  list(
    html = grep("\\.html$", all_files, value = TRUE),
    xlsx = grep("\\.xlsx$", all_files, value = TRUE),
    all = all_files
  )
}
