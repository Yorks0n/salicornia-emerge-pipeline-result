#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_file), "v2_pipeline_config.R"))

raw_dir <- file.path(V2_TEMP, "01_hog_fasta_raw")
qc_dir <- file.path(V2_TEMP, "02_hog_fasta_pre_qc")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading candidate HOGs and OrthoFinder membership...")
candidate <- read.csv(CANDIDATE_HOG_FILE, check.names = FALSE, stringsAsFactors = FALSE)
candidate_hogs <- unique(candidate$HOG[nzchar(candidate$HOG)])

n0 <- read.delim(
  N0_FILE,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)
species_columns <- setdiff(names(n0), METADATA_COLUMNS)
selected <- n0[n0$HOG %in% candidate_hogs, , drop = FALSE]

missing_hogs <- setdiff(candidate_hogs, selected$HOG)
if (length(missing_hogs) > 0L) {
  warning("Candidate HOGs absent from N0.tsv: ", paste(missing_hogs, collapse = ", "))
}

split_ids <- function(value) {
  if (is.na(value) || !nzchar(trimws(value))) return(character())
  ids <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  ids[nzchar(ids)]
}

member_parts <- list()
part_i <- 0L
for (row_i in seq_len(nrow(selected))) {
  for (species in species_columns) {
    ids <- split_ids(selected[[species]][[row_i]])
    if (length(ids) == 0L) next
    part_i <- part_i + 1L
    member_parts[[part_i]] <- data.frame(
      HOG = selected$HOG[[row_i]],
      OG = selected$OG[[row_i]],
      species = species,
      GID = ids,
      sequence_id = paste(species, ids, sep = "|"),
      stringsAsFactors = FALSE
    )
  }
}
members <- bind_rows(member_parts)

message("Loading one protein FASTA per species...")
proteomes <- setNames(vector("list", length(species_columns)), species_columns)
for (species in species_columns) {
  fasta_file <- file.path(PROTEIN_DIR, paste0(species, ".fa"))
  if (!file.exists(fasta_file)) stop("Missing protein FASTA: ", fasta_file)
  aa <- readAAStringSet(fasta_file)
  names(aa) <- sub("[[:space:]].*$", "", names(aa))
  if (anyDuplicated(names(aa))) stop("Duplicate protein IDs in ", fasta_file)
  proteomes[[species]] <- aa
}

members$found_in_fasta <- FALSE
members$sequence_length <- NA_integer_
members$x_fraction <- NA_real_
members$internal_stop <- NA

sequence_cache <- vector("list", nrow(members))
for (i in seq_len(nrow(members))) {
  aa_set <- proteomes[[members$species[[i]]]]
  gid <- members$GID[[i]]
  if (!(gid %in% names(aa_set))) next
  sequence <- as.character(aa_set[[gid]])
  sequence_cache[[i]] <- sequence
  members$found_in_fasta[[i]] <- TRUE
  members$sequence_length[[i]] <- nchar(sequence)
  aa_chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  members$x_fraction[[i]] <- if (length(aa_chars)) mean(aa_chars == "X") else NA_real_
  members$internal_stop[[i]] <- grepl("\\*", sub("\\*$", "", sequence))
}

# The HOG reference length is the median of species medians, so an expanded
# species does not dominate the length baseline.
species_length <- members %>%
  filter(found_in_fasta) %>%
  group_by(HOG, species) %>%
  summarise(species_median_length = median(sequence_length), .groups = "drop")
hog_length <- species_length %>%
  group_by(HOG) %>%
  summarise(hog_reference_length = median(species_median_length), .groups = "drop")
members <- left_join(members, hog_length, by = "HOG")
members$length_ratio <- members$sequence_length / members$hog_reference_length
members$pre_qc_pass <- with(
  members,
  found_in_fasta &
    !internal_stop &
    x_fraction <= PRE_QC_MAX_X_FRACTION &
    length_ratio >= PRE_QC_MIN_LENGTH_RATIO &
    length_ratio <= PRE_QC_MAX_LENGTH_RATIO
)
members$pre_qc_pass[is.na(members$pre_qc_pass)] <- FALSE

message("Writing raw and pre-QC HOG FASTA files...")
for (hog in unique(members$HOG)) {
  idx_raw <- which(members$HOG == hog & members$found_in_fasta)
  idx_qc <- which(members$HOG == hog & members$pre_qc_pass)

  if (length(idx_raw) > 0L) {
    raw_set <- AAStringSet(unlist(sequence_cache[idx_raw], use.names = FALSE))
    names(raw_set) <- members$sequence_id[idx_raw]
    writeXStringSet(raw_set, file.path(raw_dir, paste0(hog, ".fa")))
  }
  if (length(idx_qc) > 0L) {
    qc_set <- AAStringSet(unlist(sequence_cache[idx_qc], use.names = FALSE))
    names(qc_set) <- members$sequence_id[idx_qc]
    writeXStringSet(qc_set, file.path(qc_dir, paste0(hog, ".fa")))
  }
}

hog_summary <- members %>%
  group_by(HOG, OG) %>%
  summarise(
    n_members = n(),
    n_species = n_distinct(species),
    n_missing_fasta = sum(!found_in_fasta),
    n_pre_qc = sum(pre_qc_pass),
    n_species_pre_qc = n_distinct(species[pre_qc_pass]),
    n_target_pre_qc = sum(species == TARGET_SPECIES & pre_qc_pass),
    hog_reference_length = first(hog_reference_length),
    ready_for_alignment = n_pre_qc >= 2L & n_target_pre_qc >= 1L,
    .groups = "drop"
  )

write.csv(members, file.path(V2_TEMP, "member_metadata.csv"), row.names = FALSE, na = "")
write.csv(hog_summary, file.path(V2_TEMP, "hog_pre_qc_summary.csv"), row.names = FALSE, na = "")
write.csv(
  members[!members$found_in_fasta, c("HOG", "species", "GID")],
  file.path(V2_TEMP, "missing_protein_ids.csv"),
  row.names = FALSE,
  na = ""
)

message(
  "Prepared ", sum(hog_summary$ready_for_alignment), " / ", nrow(hog_summary),
  " HOGs for alignment. Outputs: ", V2_TEMP
)

