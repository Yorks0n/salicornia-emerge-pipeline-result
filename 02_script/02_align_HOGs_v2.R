#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_file), "v2_pipeline_config.R"))

input_dir <- file.path(V2_TEMP, "02_hog_fasta_pre_qc")
initial_dir <- file.path(V2_TEMP, "03_mafft_initial")
final_fasta_dir <- file.path(V2_TEMP, "04_hog_fasta_final")
mafft_dir <- file.path(V2_TEMP, "05_mafft_final")
trim_dir <- file.path(V2_TEMP, "06_trimal")
colmap_dir <- file.path(V2_TEMP, "07_column_maps")
log_dir <- file.path(V2_TEMP, "logs")
for (path in c(initial_dir, final_fasta_dir, mafft_dir, trim_dir, colmap_dir, log_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

metadata <- read.csv(file.path(V2_TEMP, "member_metadata.csv"), stringsAsFactors = FALSE)
summary_pre <- read.csv(file.path(V2_TEMP, "hog_pre_qc_summary.csv"), stringsAsFactors = FALSE)
hogs <- summary_pre$HOG[summary_pre$ready_for_alignment]

run_mafft <- function(input, output, log_file) {
  status <- system2(
    MAFFT_BIN,
    args = c("--auto", "--thread", as.character(MAFFT_THREADS), "--quiet", input),
    stdout = output,
    stderr = log_file
  )
  if (!identical(status, 0L) || !file.exists(output) || file.info(output)$size == 0L) {
    stop("MAFFT failed for ", basename(input), ". See ", log_file)
  }
}

parse_colmap <- function(report_file, expected_width) {
  text <- paste(readLines(report_file, warn = FALSE), collapse = " ")
  text <- sub("^.*#ColumnsMap[[:space:]]*", "", text)
  values <- as.integer(trimws(strsplit(text, ",", fixed = TRUE)[[1]]))
  values <- values[!is.na(values)]
  if (length(values) != expected_width) {
    stop(
      "trimAl column-map length mismatch for ", basename(report_file),
      ": expected ", expected_width, ", observed ", length(values)
    )
  }
  values
}

qc_records <- list()
final_members <- list()
colmap_records <- list()
q_i <- f_i <- c_i <- 0L

for (hog in hogs) {
  message("Aligning ", hog)
  input <- file.path(input_dir, paste0(hog, ".fa"))
  initial_mafft <- file.path(initial_dir, paste0(hog, ".mafft.fa"))
  initial_log <- file.path(log_dir, paste0(hog, ".mafft_initial.log"))
  run_mafft(input, initial_mafft, initial_log)

  initial_alignment <- readAAStringSet(initial_mafft)
  initial_matrix <- as.matrix(initial_alignment)
  initial_species <- sub("\\|.*$", "", rownames(initial_matrix))
  # Define the alignment core with species-equal occupancy. Terminal extensions
  # and copy-rich species therefore do not inflate the alignment denominator.
  species_occupancy <- sapply(
    unique(initial_species),
    function(species) colMeans(initial_matrix[initial_species == species, , drop = FALSE] != "-")
  )
  core_columns <- rowMeans(species_occupancy) >= 0.5
  if (!any(core_columns)) core_columns[] <- TRUE
  alignment_coverage <- rowMeans(initial_matrix[, core_columns, drop = FALSE] != "-")
  keep <- alignment_coverage >= MIN_ALIGNMENT_COVERAGE

  for (sequence_id in names(initial_alignment)) {
    q_i <- q_i + 1L
    q_idx <- match(sequence_id, names(alignment_coverage))
    qc_records[[q_i]] <- data.frame(
      HOG = hog,
      sequence_id = sequence_id,
      species_equal_core_columns = sum(core_columns),
      initial_core_coverage = unname(alignment_coverage[[q_idx]]),
      alignment_qc_pass = unname(keep[[q_idx]]),
      stringsAsFactors = FALSE
    )
  }

  retained <- initial_alignment[keep]
  retained_species <- sub("\\|.*$", "", names(retained))
  if (length(retained) < 2L || !TARGET_SPECIES %in% retained_species) {
    warning("Skipping ", hog, ": too few sequences or no target after alignment QC.")
    next
  }

  original_sequences <- readAAStringSet(input)
  final_sequences <- original_sequences[names(retained)]
  final_fasta <- file.path(final_fasta_dir, paste0(hog, ".fa"))
  writeXStringSet(final_sequences, final_fasta)

  final_mafft <- file.path(mafft_dir, paste0(hog, ".mafft.fa"))
  final_mafft_log <- file.path(log_dir, paste0(hog, ".mafft_final.log"))
  run_mafft(final_fasta, final_mafft, final_mafft_log)

  trim_file <- file.path(trim_dir, paste0(hog, ".trim.fa"))
  colmap_report <- file.path(colmap_dir, paste0(hog, ".colmap.raw.txt"))
  trim_log <- file.path(log_dir, paste0(hog, ".trimal.log"))
  status <- system2(
    TRIMAL_BIN,
    args = c("-in", final_mafft, "-out", trim_file, "-gappyout", "-colnumbering"),
    stdout = colmap_report,
    stderr = trim_log
  )
  if (!identical(status, 0L) || !file.exists(trim_file) || file.info(trim_file)$size == 0L) {
    stop("trimAl failed for ", hog, ". See ", trim_log)
  }

  trimmed <- readAAStringSet(trim_file)
  mafft_columns_0based <- parse_colmap(colmap_report, width(trimmed)[[1]])
  colmap <- data.frame(
    HOG = hog,
    trim_col_1based = seq_along(mafft_columns_0based),
    mafft_col_0based = mafft_columns_0based,
    mafft_col_1based = mafft_columns_0based + 1L
  )
  write.csv(colmap, file.path(colmap_dir, paste0(hog, ".colmap.csv")), row.names = FALSE)
  c_i <- c_i + 1L
  colmap_records[[c_i]] <- colmap

  post_matrix <- as.matrix(trimmed)
  post_coverage <- rowMeans(post_matrix != "-")
  for (sequence_id in names(final_sequences)) {
    meta_row <- metadata[metadata$sequence_id == sequence_id & metadata$HOG == hog, , drop = FALSE]
    f_i <- f_i + 1L
    final_members[[f_i]] <- data.frame(
      meta_row,
      trimmed_alignment_coverage = unname(post_coverage[[sequence_id]]),
      stringsAsFactors = FALSE
    )
  }
}

qc_df <- bind_rows(qc_records)
final_df <- bind_rows(final_members)
all_colmaps <- bind_rows(colmap_records)
write.csv(qc_df, file.path(V2_TEMP, "alignment_qc_metrics.csv"), row.names = FALSE, na = "")
write.csv(final_df, file.path(V2_TEMP, "final_member_metadata.csv"), row.names = FALSE, na = "")
write.csv(all_colmaps, file.path(V2_TEMP, "all_column_maps.csv"), row.names = FALSE)

message("Finished final alignment and trimAl mapping for ", length(unique(final_df$HOG)), " HOGs.")
