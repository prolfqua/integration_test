#!/usr/bin/env Rscript
# Compare prolfquapp DEA outputs across Docker image versions / models.
#
# Reads each run's DE_*.xlsx `diff_exp_analysis` sheet (plus
# `diff_exp_analysis_wide` for FDR), joins protein x contrast across runs, and
# reports Pearson correlation + max absolute difference for the fold-change
# (diff), p.value and FDR estimates.
#
# Usage:
#   Rscript compare_docker_versions.R [--min-cor=0.95] [--pairs=A:B,C:D] \
#       <labelA>=<DE_A.xlsx> <labelB>=<DE_B.xlsx> ...
#
# Each run is "<label>=<path to DE_*.xlsx>". Labels may not contain "=".
# --pairs  explicit "left:right" comparisons (comma-separated). If omitted,
#          consecutive runs plus first-vs-last are compared.
# --min-cor  fold-change/p.value/FDR correlations below this fail the run
#            (non-zero exit). Default 0.95.

suppressMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(min_cor = 0.95, pairs = NULL)
runs <- character()
for (a in args) {
  if (startsWith(a, "--min-cor=")) {
    opt$min_cor <- as.numeric(sub("--min-cor=", "", a))
  } else if (startsWith(a, "--pairs=")) {
    opt$pairs <- strsplit(sub("--pairs=", "", a), ",", fixed = TRUE)[[1]]
  } else {
    runs <- c(runs, a)
  }
}
if (length(runs) < 2) stop("need at least two <label>=<DE.xlsx> runs")

specs <- lapply(runs, function(a) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  list(label = kv[1], path = paste(kv[-1], collapse = "="))
})

# Long-format DEA result with FDR merged in from the wide sheet.
read_dea <- function(path) {
  long <- readxl::read_excel(path, sheet = "diff_exp_analysis")
  long <- long[, intersect(
    c("protein_Id", "contrast", "diff", "p.value", "statistic"),
    colnames(long)
  )]
  wide <- tryCatch(
    readxl::read_excel(path, sheet = "diff_exp_analysis_wide"),
    error = function(e) NULL
  )
  fdr_cols <- if (is.null(wide)) character() else grep("^FDR\\.", colnames(wide), value = TRUE)
  if (length(fdr_cols) > 0) {
    fdr_long <- wide |>
      dplyr::select(dplyr::all_of(c("protein_Id", fdr_cols))) |>
      tidyr::pivot_longer(dplyr::all_of(fdr_cols),
                          names_to = "contrast", values_to = "FDR") |>
      dplyr::mutate(contrast = sub("^FDR\\.", "", contrast))
    long <- dplyr::left_join(long, fdr_long, by = c("protein_Id", "contrast"))
  }
  long
}

dat <- lapply(specs, function(s) read_dea(s$path))
names(dat) <- vapply(specs, function(s) s$label, "")

cat("Loaded runs:\n")
for (nm in names(dat)) {
  cat(sprintf("  %-18s %d rows, contrasts: %s\n", nm, nrow(dat[[nm]]),
              paste(unique(dat[[nm]]$contrast), collapse = ", ")))
}

corr <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(a[ok]) == 0 && sd(b[ok]) == 0) return(1)
  cor(a[ok], b[ok])
}

labels <- names(dat)
if (is.null(opt$pairs)) {
  pairs <- Map(c, labels[-length(labels)], labels[-1])
  if (length(labels) > 2) pairs <- c(pairs, list(c(labels[1], labels[length(labels)])))
} else {
  pairs <- lapply(opt$pairs, function(p) strsplit(p, ":", fixed = TRUE)[[1]])
}

metrics <- c("diff", "p.value", "FDR")
cat(sprintf("\n%-34s %-12s %-7s %12s %12s %6s\n",
            "pair", "contrast", "metric", "r", "max|delta|", "n"))
cat(strrep("-", 90), "\n")
fail <- FALSE
for (p in pairs) {
  if (!all(p %in% labels)) {
    cat(sprintf("SKIP %s vs %s (unknown label)\n", p[1], p[2])); next
  }
  m <- dplyr::inner_join(dat[[p[1]]], dat[[p[2]]],
                         by = c("protein_Id", "contrast"), suffix = c(".x", ".y"))
  for (ct in sort(unique(m$contrast))) {
    mc <- m[m$contrast == ct, ]
    for (met in metrics) {
      cx <- paste0(met, ".x"); cy <- paste0(met, ".y")
      if (!all(c(cx, cy) %in% colnames(mc))) next
      r <- corr(mc[[cx]], mc[[cy]])
      md <- suppressWarnings(max(abs(mc[[cx]] - mc[[cy]]), na.rm = TRUE))
      n <- sum(is.finite(mc[[cx]]) & is.finite(mc[[cy]]))
      flag <- if (!is.na(r) && r < opt$min_cor) { fail <- TRUE; " FAIL" } else ""
      cat(sprintf("%-34s %-12s %-7s %12.6f %12.4g %6d%s\n",
                  paste0(p[1], " vs ", p[2]), ct, met, r, md, n, flag))
    }
  }
}
cat(sprintf("\nmin-cor threshold: %.3f -> %s\n", opt$min_cor,
            if (fail) "FAIL" else "PASS"))
if (fail) quit(status = 1)
