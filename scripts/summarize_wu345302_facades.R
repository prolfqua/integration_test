args <- commandArgs(trailingOnly = TRUE)
out_root <- if (length(args) >= 1) args[[1]] else file.path("test-outputs", "wu345302_facades")

status_path <- file.path(out_root, "status.tsv")
if (!file.exists(status_path)) {
  stop("Missing status file: ", status_path)
}

status <- utils::read.delim(status_path, stringsAsFactors = FALSE)
ok <- status[status$run_status == "ok", , drop = FALSE]

find_xlsx <- function(output_dir) {
  files <- list.files(
    output_dir,
    pattern = "^DE_WU345302[.]xlsx$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(NA_character_)
  }
  files[[1]]
}

find_index <- function(output_dir) {
  files <- list.files(
    output_dir,
    pattern = "^index[.]html$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(NA_character_)
  }
  files[[1]]
}

read_result <- function(model, xlsx) {
  d <- readxl::read_xlsx(xlsx, sheet = "diff_exp_analysis")
  fdr <- if ("FDR" %in% names(d)) {
    d$FDR
  } else if ("BFDR" %in% names(d)) {
    d$BFDR
  } else if ("adj.P.Val" %in% names(d)) {
    d$adj.P.Val
  } else {
    rep(NA_real_, nrow(d))
  }
  diff <- if ("diff" %in% names(d)) {
    d$diff
  } else if ("log2_EFCs" %in% names(d)) {
    d$log2_EFCs
  } else if ("logFC" %in% names(d)) {
    d$logFC
  } else {
    rep(NA_real_, nrow(d))
  }
  p_value <- if ("p.value" %in% names(d)) {
    d$p.value
  } else if ("P.Value" %in% names(d)) {
    d$P.Value
  } else {
    rep(NA_real_, nrow(d))
  }
  model_name <- if ("modelName" %in% names(d)) {
    paste(unique(d$modelName), collapse = "|")
  } else {
    NA_character_
  }

  list(
    summary = data.frame(
      model = model,
      rows = nrow(d),
      finite_FDR = sum(is.finite(fdr)),
      FDR_lt_0_05 = sum(fdr < 0.05, na.rm = TRUE),
      finite_diff = sum(is.finite(diff)),
      modelName = model_name
    ),
    data = data.frame(
      protein_Id = d$protein_Id,
      model = model,
      FDR = fdr,
      diff = diff,
      p.value = p_value
    )
  )
}

results <- list()
summaries <- list()
for (i in seq_len(nrow(ok))) {
  model <- ok$model[[i]]
  xlsx <- find_xlsx(ok$output_dir[[i]])
  if (is.na(xlsx)) {
    next
  }
  res <- read_result(model, xlsx)
  summaries[[model]] <- res$summary
  results[[model]] <- res$data
}

model_summary <- do.call(rbind, summaries)
model_summary$index_html <- vapply(
  model_summary$model,
  function(model) {
    find_index(status$output_dir[match(model, status$model)])
  },
  character(1)
)
utils::write.table(
  model_summary,
  file.path(out_root, "model_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

wide_metric <- function(metric) {
  Reduce(
    function(x, y) merge(x, y, by = "protein_Id", all = TRUE),
    lapply(names(results), function(model) {
      z <- results[[model]][, c("protein_Id", metric)]
      names(z)[[2]] <- model
      z
    })
  )
}

pairwise <- data.frame()
ref <- "limma_impute"
if (ref %in% names(results)) {
  fdr_wide <- wide_metric("FDR")
  diff_wide <- wide_metric("diff")
  for (model in setdiff(names(results), ref)) {
    keep_fdr <- is.finite(fdr_wide[[ref]]) & is.finite(fdr_wide[[model]])
    keep_diff <- is.finite(diff_wide[[ref]]) & is.finite(diff_wide[[model]])
    pairwise <- rbind(
      pairwise,
      data.frame(
        reference = ref,
        model = model,
        n_FDR = sum(keep_fdr),
        pearson_FDR = cor(fdr_wide[[ref]][keep_fdr], fdr_wide[[model]][keep_fdr]),
        spearman_FDR = cor(
          fdr_wide[[ref]][keep_fdr],
          fdr_wide[[model]][keep_fdr],
          method = "spearman"
        ),
        pearson_neglog10_FDR = cor(
          -log10(pmax(fdr_wide[[ref]][keep_fdr], .Machine$double.xmin)),
          -log10(pmax(fdr_wide[[model]][keep_fdr], .Machine$double.xmin))
        ),
        n_diff = sum(keep_diff),
        pearson_diff = cor(diff_wide[[ref]][keep_diff], diff_wide[[model]][keep_diff])
      )
    )
  }
}
utils::write.table(
  pairwise,
  file.path(out_root, "pairwise_vs_limma_impute.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

failures <- status[status$run_status != "ok", , drop = FALSE]
if (nrow(failures) > 0) {
  failure_summary <- lapply(seq_len(nrow(failures)), function(i) {
    model <- failures$model[[i]]
    log_file <- file.path(out_root, "logs", paste0(model, ".stdout.log"))
    lines <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character()
    error_lines <- grep("ERROR|Error|Execution halted", lines, value = TRUE)
    data.frame(
      model = model,
      exit_code = failures$exit_code[[i]],
      log_file = log_file,
      error = paste(tail(error_lines, 5), collapse = " | ")
    )
  })
  failure_summary <- do.call(rbind, failure_summary)
  utils::write.table(
    failure_summary,
    file.path(out_root, "failures.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

message("Wrote: ", file.path(out_root, "model_summary.tsv"))
message("Wrote: ", file.path(out_root, "pairwise_vs_limma_impute.tsv"))
if (nrow(failures) > 0) {
  message("Wrote: ", file.path(out_root, "failures.tsv"))
}
