#!/usr/bin/env Rscript
# Compare prolfquapp DEA outputs across Docker image versions.
#
# Each version is run (outside this script) on the SAME fixture with matched
# processing options (model=prolfqua, aggregate=medpolish, transform=robscale,
# FDR=0.1). This script reads the per-version DE_*.xlsx `diff_exp_analysis`
# sheet, joins protein x contrast across versions, and reports Pearson
# correlation and max absolute difference for the fold-change (diff),
# p.value and FDR estimates.
#
# Usage:
#   Rscript scripts/compare_docker_versions.R <label1>=<DE1.xlsx> <label2>=<DE2.xlsx> ...
# Versions are compared pairwise in the order given (v1-v2, v2-v3, v1-v3 ...).

suppressMessages({
  library(readxl)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("need at least two <label>=<DE.xlsx> arguments")
}

parse_arg <- function(a) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  list(label = kv[1], path = paste(kv[-1], collapse = "="))
}
specs <- lapply(args, parse_arg)

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
  if (!is.null(wide)) {
    fdr_cols <- grep("^FDR\\.", colnames(wide), value = TRUE)
    if (length(fdr_cols) > 0) {
      fdr_long <- wide |>
        dplyr::select(dplyr::all_of(c("protein_Id", fdr_cols))) |>
        tidyr::pivot_longer(
          dplyr::all_of(fdr_cols),
          names_to = "contrast", values_to = "FDR"
        ) |>
        dplyr::mutate(contrast = sub("^FDR\\.", "", contrast))
      long <- dplyr::left_join(long, fdr_long, by = c("protein_Id", "contrast"))
    }
  }
  long
}

dat <- lapply(specs, function(s) read_dea(s$path))
names(dat) <- vapply(specs, function(s) s$label, "")

cat("Loaded versions:\n")
for (nm in names(dat)) {
  cat(sprintf("  %-8s %d rows, contrasts: %s\n", nm, nrow(dat[[nm]]),
              paste(unique(dat[[nm]]$contrast), collapse = ", ")))
}

metrics <- c("diff", "p.value", "FDR")
corr <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(a[ok]) == 0 && sd(b[ok]) == 0) return(1)
  cor(a[ok], b[ok])
}

# Pairwise comparisons: consecutive pairs + first-vs-last.
labels <- names(dat)
pairs <- list()
for (i in seq_len(length(labels) - 1)) {
  pairs[[length(pairs) + 1]] <- c(labels[i], labels[i + 1])
}
if (length(labels) > 2) {
  pairs[[length(pairs) + 1]] <- c(labels[1], labels[length(labels)])
}

cat(sprintf("\n%-22s %-12s %-7s %12s %12s %6s\n",
            "pair", "contrast", "metric", "r", "max|delta|", "n"))
cat(strrep("-", 78), "\n")
for (p in pairs) {
  a <- dat[[p[1]]]; b <- dat[[p[2]]]
  m <- dplyr::inner_join(a, b, by = c("protein_Id", "contrast"),
                         suffix = c(".x", ".y"))
  for (ct in sort(unique(m$contrast))) {
    mc <- m[m$contrast == ct, ]
    for (met in metrics) {
      cx <- paste0(met, ".x"); cy <- paste0(met, ".y")
      if (!all(c(cx, cy) %in% colnames(mc))) next
      r <- corr(mc[[cx]], mc[[cy]])
      md <- suppressWarnings(max(abs(mc[[cx]] - mc[[cy]]), na.rm = TRUE))
      cat(sprintf("%-22s %-12s %-7s %12.6f %12.4g %6d\n",
                  paste0(p[1], " vs ", p[2]), ct, met, r, md,
                  sum(is.finite(mc[[cx]]) & is.finite(mc[[cy]]))))
    }
  }
}
