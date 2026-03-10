#!/usr/bin/env Rscript
# create_test_fixtures.R
#
# One-time script to generate small fixture datasets for integration testing.
# Each fixture contains a subset of ~100 proteins from real data, in the
# native format each prolfquapp preprocessor expects.
#
# Usage (from integration_test/):
#   Rscript R/create_test_fixtures.R
#
# Prerequisites:
#   - prolfquadata package installed (IonStar MaxQuant + MSFragger ZIPs)
#   - seqinr package for FASTA subsetting
#   - Internet connection (downloads PTM data from Zenodo on first run)

library(seqinr)

FIXTURE_DIR <- "fixtures"
CACHE_DIR <- ".cache"
N_PROTEINS <- 100
set.seed(42)

ZENODO_URL <- "https://zenodo.org/records/15879865/files/PTM_experiment_FP_22_Maculins_and_QC.zip"
ZENODO_MD5 <- "b357973f99f2dc743a5984480ec351b4"


# =============================================================================
# Zenodo download + cache
# =============================================================================
get_ptm_data_dir <- function() {
  ptm_dir <- file.path(CACHE_DIR, "PTM_example")
  if (dir.exists(ptm_dir)) {
    message("Using cached PTM data: ", ptm_dir)
    return(ptm_dir)
  }

  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  zip_path <- file.path(CACHE_DIR, basename(ZENODO_URL))

  message("Downloading PTM data from Zenodo (~91 MB)...")
  download.file(ZENODO_URL, zip_path, mode = "wb")

  # Verify checksum
  actual_md5 <- tools::md5sum(zip_path)
  if (actual_md5 != ZENODO_MD5) {
    unlink(zip_path)
    stop("MD5 mismatch! Expected ", ZENODO_MD5, " got ", actual_md5)
  }
  message("  MD5 OK")

  message("Extracting...")
  unzip(zip_path, exdir = CACHE_DIR)
  unlink(zip_path)

  # The ZIP extracts to PTM_example/
  if (!dir.exists(ptm_dir)) {
    # Check what was actually extracted
    extracted <- list.dirs(CACHE_DIR, recursive = FALSE)
    stop("Expected PTM_example/ after extraction, found: ",
         paste(basename(extracted), collapse = ", "))
  }

  message("  PTM data cached at: ", ptm_dir)
  return(ptm_dir)
}


# =============================================================================
# Generate dataset_with_contrasts.tsv using prolfqua
# =============================================================================
generate_dataset_with_contrasts <- function(output_path) {
  # 22 TMT channels: 2 genotypes x 3 timepoints x 3-4 replicates
  channels <- c(
    "KO_Early_1", "KO_Early_2", "KO_Early_3", "KO_Early_4",
    "KO_Late_1", "KO_Late_2", "KO_Late_3", "KO_Late_4",
    "KO_Uninfect_1", "KO_Uninfect_2", "KO_Uninfect_3",
    "WT_Early_1", "WT_Early_2", "WT_Early_3", "WT_Early_4",
    "WT_Late_1", "WT_Late_2", "WT_Late_3", "WT_Late_4",
    "WT_Uninfect_1", "WT_Uninfect_3", "WT_Uninfect_4"
  )
  annot <- data.frame(channel = channels, Name = channels,
                      stringsAsFactors = FALSE)
  annot <- tidyr::separate(annot, "Name", c("Genotype", "Timepoint", NA),
                           sep = "_", remove = FALSE)
  annot2 <- prolfqua::annotation_add_contrasts(
    annot, primary_col = "Genotype", secondary_col = "Timepoint",
    decreasing = TRUE, interactions = FALSE
  )$annot
  readr::write_tsv(annot2, output_path)
  message("  Generated: ", output_path)
}


# =============================================================================
# Helper: subset FASTA to matching protein IDs
# =============================================================================
subset_fasta <- function(fasta_path, protein_ids, output_path) {
  seqs <- seqinr::read.fasta(fasta_path, seqtype = "AA", whole.header = TRUE)
  seq_names <- names(seqs)

  get_accession <- function(name) {
    if (grepl("^(sp|tr)\\|", name)) {
      sub("^(sp|tr)\\|([^|]+)\\|.*", "\\2", name)
    } else {
      sub("\\s.*", "", name)
    }
  }
  accessions <- sapply(seq_names, get_accession, USE.NAMES = FALSE)

  keep <- seq_names %in% protein_ids |
    accessions %in% protein_ids |
    sapply(seq_names, function(n) any(grepl(n, protein_ids, fixed = TRUE)))

  if (sum(keep) == 0) {
    keep <- sapply(seq_names, function(n) {
      any(sapply(protein_ids, function(pid) grepl(pid, n, fixed = TRUE)))
    })
  }

  message("  FASTA: keeping ", sum(keep), " of ", length(seqs), " sequences")
  if (sum(keep) > 0) {
    seqinr::write.fasta(
      sequences = seqs[keep],
      names = sapply(seqs[keep], attr, "Annot"),
      file.out = output_path,
      as.string = TRUE
    )
  }
}


# =============================================================================
# 1a. MaxQuant IonStar fixture
# =============================================================================
create_maxquant_fixture <- function() {
  message("=== Creating MaxQuant IonStar fixture ===")
  outdir <- file.path(FIXTURE_DIR, "maxquant_ionstar")
  dir.create(file.path(outdir, "txt"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "fasta"), showWarnings = FALSE)

  annot_file <- system.file("quantdata/annotation_Ionstar2018_PXD003881.xlsx",
                            package = "prolfquadata")
  annot <- readxl::read_xlsx(annot_file)
  message("  Annotation: ", nrow(annot), " samples")

  zip_file <- system.file("quantdata/MAXQuant_IonStar2018_PXD003881.zip",
                          package = "prolfquadata")
  tmp <- tempdir()
  unzip(zip_file, files = "peptides.txt", exdir = tmp)
  peptides <- read.csv(file.path(tmp, "peptides.txt"), sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
  message("  peptides.txt: ", nrow(peptides), " rows, ",
          length(unique(peptides$`Leading razor protein`)), " proteins")

  prot_counts <- sort(table(peptides$`Leading razor protein`), decreasing = TRUE)
  prot_names <- names(prot_counts)
  prot_names <- prot_names[!grepl("^REV_|^CON_|^zz|^rev_", prot_names)]
  selected_prots <- head(prot_names, N_PROTEINS)

  peptides_sub <- peptides[peptides$`Leading razor protein` %in% selected_prots, ]
  message("  Subset: ", nrow(peptides_sub), " peptides for ", length(selected_prots), " proteins")
  write.table(peptides_sub, file.path(outdir, "txt", "peptides.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  ds <- data.frame(
    raw.file = annot$raw.file,
    Name = paste0(annot$sample, "_", annot$run_ID),
    Group = annot$sample,
    Subject = annot$run_ID,
    CONTROL = ifelse(annot$sample == "b", "C",
              ifelse(annot$sample == "e", "T", NA)),
    stringsAsFactors = FALSE
  )
  ds <- ds[!is.na(ds$CONTROL), ]
  message("  Annotation: ", nrow(ds), " samples (groups b vs e)")
  write.csv(ds, file.path(outdir, "dataset.csv"), row.names = FALSE)

  fasta_file <- system.file("fastaDBs/uniprot-proteome_UP000005640_reviewed_yes.fasta.gz",
                            package = "prolfquadata")
  if (file.exists(fasta_file)) {
    subset_fasta(fasta_file, selected_prots, file.path(outdir, "fasta", "test.fasta"))
  } else {
    message("  WARNING: No FASTA found in prolfquadata, creating minimal dummy")
    lines <- sapply(selected_prots, function(p) {
      paste0(">", p, " Dummy protein\n", paste(rep("A", 100), collapse = ""))
    })
    writeLines(unlist(lines), file.path(outdir, "fasta", "test.fasta"))
  }

  config <- list(
    path = ".",
    zipdir_name = "DEA_test_maxquant",
    prefix = "DEA",
    software = "prolfquapp.MAXQUANT",
    project_spec = list(
      input_URL = "", workunit_Id = "test_maxquant",
      order_Id = "", project_name = "integration_test", project_Id = ""
    ),
    processing_options = list(
      model = "prolfqua", model_missing = TRUE, interaction = FALSE,
      nr_peptides = 1,
      pattern_contaminants = "^zz|^CON|Cont_",
      pattern_decoys = "^REV_|^rev_",
      remove_decoys = FALSE, remove_cont = FALSE,
      FDR_threshold = 0.1, diff_threshold = 1.0,
      aggregate = "medpolish", transform = "robscale"
    ),
    ext_reader = list(dataset = list(), extra_args = "list()",
                      preprocess = list(), get_files = list()),
    group = "G_"
  )
  yaml::write_yaml(config, file.path(outdir, "config.yaml"))
  message("  Done: ", outdir)
}


# =============================================================================
# 1b. FragPipe MSstats fixture
# =============================================================================
create_msstats_fixture <- function() {
  message("\n=== Creating FragPipe MSstats IonStar fixture ===")
  outdir <- file.path(FIXTURE_DIR, "fragpipe_ionstar")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "fasta"), showWarnings = FALSE)

  annot_file <- system.file("quantdata/annotation_Ionstar2018_PXD003881.xlsx",
                            package = "prolfquadata")
  annot <- readxl::read_xlsx(annot_file)

  zip_file <- system.file("quantdata/MSFragger_IonStar2018_PXD003881.zip",
                          package = "prolfquadata")
  tmp <- tempdir()
  unzip(zip_file, files = "IonstarWithMSFragger/MSstats.csv", exdir = tmp)
  msstats <- readr::read_csv(file.path(tmp, "IonstarWithMSFragger", "MSstats.csv"),
                             show_col_types = FALSE)
  message("  MSstats.csv: ", nrow(msstats), " rows, ",
          length(unique(msstats$ProteinName)), " proteins")

  prot_counts <- sort(table(msstats$ProteinName), decreasing = TRUE)
  prot_names <- names(prot_counts)
  prot_names <- prot_names[!grepl("^REV_|^CON_|^zz|^rev_", prot_names)]
  selected_prots <- head(prot_names, 50)

  msstats_sub <- msstats[msstats$ProteinName %in% selected_prots, ]
  message("  Subset: ", nrow(msstats_sub), " rows for ", length(selected_prots), " proteins")
  readr::write_csv(msstats_sub, file.path(outdir, "MSstats.csv"))

  runs <- sort(unique(msstats_sub$Run))
  ds <- data.frame(raw.file = runs, stringsAsFactors = FALSE)
  ds$sample <- toupper(sub(".*_human_ecoli_([A-Za-z])_.*", "\\1", ds$raw.file))
  ds$run_num <- sub("^B03_(\\d+)_.*", "\\1", ds$raw.file)
  ds$Name <- paste0(ds$sample, "_", ds$run_num)
  ds$Group <- ds$sample
  ds$Subject <- ds$run_num
  ds$CONTROL <- ifelse(ds$sample == "B", "C", ifelse(ds$sample == "E", "T", NA))
  ds <- ds[!is.na(ds$CONTROL), ]
  ds <- ds[, c("raw.file", "Name", "Group", "Subject", "CONTROL")]
  message("  Annotation: ", nrow(ds), " samples (groups B vs E)")
  write.csv(ds, file.path(outdir, "dataset.csv"), row.names = FALSE)

  fasta_file <- system.file("fastaDBs/uniprot-proteome_UP000005640_reviewed_yes.fasta.gz",
                            package = "prolfquadata")
  if (file.exists(fasta_file)) {
    subset_fasta(fasta_file, selected_prots, file.path(outdir, "fasta", "test.fasta"))
  } else {
    lines <- sapply(selected_prots, function(p) {
      paste0(">", p, " Dummy protein\n", paste(rep("A", 100), collapse = ""))
    })
    writeLines(unlist(lines), file.path(outdir, "fasta", "test.fasta"))
  }

  config <- list(
    path = ".",
    zipdir_name = "DEA_test_msstats",
    prefix = "DEA",
    software = "prolfquapp.MSSTATS",
    project_spec = list(
      input_URL = "", workunit_Id = "test_msstats",
      order_Id = "", project_name = "integration_test", project_Id = ""
    ),
    processing_options = list(
      model = "prolfqua", model_missing = TRUE, interaction = FALSE,
      nr_peptides = 1,
      pattern_contaminants = "^zz|^CON|Cont_",
      pattern_decoys = "^REV_|^rev_",
      remove_decoys = FALSE, remove_cont = FALSE,
      FDR_threshold = 0.1, diff_threshold = 1.0,
      aggregate = "medpolish", transform = "robscale"
    ),
    ext_reader = list(dataset = list(), extra_args = "list()",
                      preprocess = list(), get_files = list()),
    group = "G_"
  )
  yaml::write_yaml(config, file.path(outdir, "config.yaml"))
  message("  Done: ", outdir)
}


# =============================================================================
# 1c. FragPipe TMT total proteome fixture
# =============================================================================
create_fp_tmt_fixture <- function(ptm_dir) {
  message("\n=== Creating FragPipe TMT total proteome fixture ===")
  outdir <- file.path(FIXTURE_DIR, "fp_tmt_total")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "fasta"), showWarnings = FALSE)

  # FP_TMT preprocessor finds psm.tsv files recursively, so preserve
  # the plex subdirectory structure (p1/psm.tsv, p2/psm.tsv)
  psm_files <- dir(file.path(ptm_dir, "data_total", "FP_22"),
                   pattern = "psm\\.tsv$", recursive = TRUE, full.names = TRUE)
  stopifnot(length(psm_files) > 0)

  # Select proteins using the first plex, then subset all plexes
  psm1 <- readr::read_tsv(psm_files[1], show_col_types = FALSE)
  prot_counts <- sort(table(psm1$Protein), decreasing = TRUE)
  prot_names <- names(prot_counts)
  prot_names <- prot_names[!grepl("^REV_|^CON_|^zz|^rev_", prot_names)]
  selected_prots <- head(prot_names, N_PROTEINS)

  total_psms <- 0
  for (psm_file in psm_files) {
    # Preserve subdirectory: e.g. p1/psm.tsv -> outdir/p1/psm.tsv
    rel_path <- sub(paste0(normalizePath(file.path(ptm_dir, "data_total", "FP_22")), "/"), "",
                    normalizePath(psm_file))
    out_path <- file.path(outdir, rel_path)
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

    psm <- readr::read_tsv(psm_file, show_col_types = FALSE)
    psm_sub <- psm[psm$Protein %in% selected_prots, ]
    total_psms <- total_psms + nrow(psm_sub)
    readr::write_tsv(psm_sub, out_path)
    message("  ", rel_path, ": ", nrow(psm_sub), " PSMs")
  }
  message("  Total: ", total_psms, " PSMs for ", length(selected_prots), " proteins")

  # Generate dataset_with_contrasts.tsv
  generate_dataset_with_contrasts(file.path(outdir, "dataset_with_contrasts.tsv"))

  # Subset FASTA
  fasta_src <- file.path(ptm_dir, "data_total", "FP_22",
                         "p37688_db3_MusNShigella_20250219.fasta")
  if (file.exists(fasta_src)) {
    accessions <- sub("^(sp|tr)\\|([^|]+)\\|.*", "\\2", selected_prots)
    subset_fasta(fasta_src, c(selected_prots, accessions),
                 file.path(outdir, "fasta", "test.fasta"))
  }

  config <- list(
    path = ".",
    zipdir_name = "DEA_test_fp_tmt",
    prefix = "DEA",
    software = "prolfquapp.FP_TMT",
    project_spec = list(
      input_URL = "", workunit_Id = "test_fp_tmt",
      order_Id = "", project_name = "integration_test", project_Id = ""
    ),
    processing_options = list(
      model = "prolfqua", model_missing = TRUE, interaction = FALSE,
      nr_peptides = 1,
      pattern_contaminants = "^zz|^CON|Cont_",
      pattern_decoys = "^REV_|^rev_",
      remove_decoys = FALSE, remove_cont = FALSE,
      FDR_threshold = 0.1, diff_threshold = 1.0,
      aggregate = "medpolish", transform = "vsn"
    ),
    ext_reader = list(dataset = list(), extra_args = "list()",
                      preprocess = list(), get_files = list()),
    group = "G_"
  )
  yaml::write_yaml(config, file.path(outdir, "config.yaml"))
  message("  Done: ", outdir)
}


# =============================================================================
# 1d. FragPipe singlesite phospho fixture
# =============================================================================
create_fp_singlesite_fixture <- function(ptm_dir) {
  message("\n=== Creating FragPipe singlesite phospho fixture ===")
  outdir <- file.path(FIXTURE_DIR, "fp_singlesite_phospho")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "fasta"), showWarnings = FALSE)

  site_file <- file.path(ptm_dir, "data_ptm", "FP_22",
                         "abundance_single-site_None.tsv")
  stopifnot(file.exists(site_file))
  site_data <- readr::read_tsv(site_file, show_col_types = FALSE)
  message("  abundance_single-site_None.tsv: ", nrow(site_data), " rows, ",
          length(unique(site_data$ProteinID)), " proteins")

  prot_counts <- sort(table(site_data$ProteinID), decreasing = TRUE)
  prot_names <- names(prot_counts)
  prot_names <- prot_names[!grepl("^REV_|^CON_|^zz|^rev_", prot_names)]
  selected_prots <- head(prot_names, 50)

  site_sub <- site_data[site_data$ProteinID %in% selected_prots, ]
  message("  Subset: ", nrow(site_sub), " sites for ", length(selected_prots), " proteins")
  readr::write_tsv(site_sub, file.path(outdir, "abundance_single-site_None.tsv"))

  # Same annotation as total proteome (same 22 samples, same contrasts)
  generate_dataset_with_contrasts(file.path(outdir, "dataset_with_contrasts.tsv"))

  # Subset FASTA
  fasta_src <- file.path(ptm_dir, "data_ptm", "FP_22",
                         "p37688_db3_MusNShigella_20250219.fasta")
  if (file.exists(fasta_src)) {
    subset_fasta(fasta_src, selected_prots,
                 file.path(outdir, "fasta", "test.fasta"))
  }

  config <- list(
    path = ".",
    zipdir_name = "DEA_test_fp_singlesite",
    prefix = "DEA",
    software = "prolfquappPTMreaders.FP_singlesite",
    project_spec = list(
      input_URL = "", workunit_Id = "test_fp_singlesite",
      order_Id = "", project_name = "integration_test", project_Id = ""
    ),
    processing_options = list(
      model = "prolfqua", model_missing = TRUE, interaction = FALSE,
      nr_peptides = 1,
      pattern_contaminants = "^zz|^CON|Cont_",
      pattern_decoys = "^REV_|^rev_",
      remove_decoys = FALSE, remove_cont = FALSE,
      FDR_threshold = 0.1, diff_threshold = 1.0,
      aggregate = "medpolish", transform = "vsn"
    ),
    ext_reader = list(dataset = list(), extra_args = "list()",
                      preprocess = list(), get_files = list()),
    group = "G_"
  )
  yaml::write_yaml(config, file.path(outdir, "config.yaml"))
  message("  Done: ", outdir)
}


# =============================================================================
# Main
# =============================================================================
message("Creating integration test fixtures in: ", FIXTURE_DIR)
message("Working directory: ", getwd())
dir.create(FIXTURE_DIR, recursive = TRUE, showWarnings = FALSE)

create_maxquant_fixture()
create_msstats_fixture()

# PTM fixtures: download from Zenodo (cached in .cache/)
ptm_dir <- get_ptm_data_dir()
create_fp_tmt_fixture(ptm_dir)
create_fp_singlesite_fixture(ptm_dir)

message("\n=== All fixtures created successfully ===")
message("Fixture directory: ", normalizePath(FIXTURE_DIR))
