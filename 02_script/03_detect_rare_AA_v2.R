#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_file), "v2_pipeline_config.R"))
if (SLIDING_WINDOW_AA %% 2L != 1L) stop("SLIDING_WINDOW_AA must be an odd integer.")
SLIDING_FLANK <- as.integer((SLIDING_WINDOW_AA - 1L) / 2L)

final_fasta_dir <- file.path(V2_TEMP, "04_hog_fasta_final")
mafft_dir <- file.path(V2_TEMP, "05_mafft_final")
trim_dir <- file.path(V2_TEMP, "06_trimal")
colmap_dir <- file.path(V2_TEMP, "07_column_maps")
site_dir <- file.path(V2_TEMP, "08_site_evidence")
dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(file.path(V2_TEMP, "final_member_metadata.csv"), stringsAsFactors = FALSE)
hogs <- unique(metadata$HOG)

is_valid_aa <- function(x) x %in% STANDARD_AA

compute_background_metrics <- function(alignment_matrix, sequence_species, background_species) {
  n_columns <- ncol(alignment_matrix)
  n_background <- length(background_species)
  if (n_background == 0L) stop("No background species are available.")

  species_non_gap <- sapply(
    background_species,
    function(species) {
      colMeans(alignment_matrix[sequence_species == species, , drop = FALSE] != "-")
    }
  )
  if (is.null(dim(species_non_gap))) species_non_gap <- matrix(species_non_gap, ncol = 1L)

  consensus_aa <- rep(NA_character_, n_columns)
  consensus_frequency <- rep(NA_real_, n_columns)
  species_valid <- integer(n_columns)

  for (column in seq_len(n_columns)) {
    species_frequencies <- list()
    for (species in background_species) {
      residues <- alignment_matrix[sequence_species == species, column]
      residues <- residues[is_valid_aa(residues)]
      if (length(residues) == 0L) next
      species_frequencies[[species]] <- table(factor(residues, levels = STANDARD_AA)) / length(residues)
    }
    species_valid[[column]] <- length(species_frequencies)
    if (length(species_frequencies) == 0L) next
    weighted_frequency <- Reduce("+", species_frequencies) / length(species_frequencies)
    consensus_index <- which.max(weighted_frequency)
    consensus_aa[[column]] <- names(weighted_frequency)[[consensus_index]]
    consensus_frequency[[column]] <- unname(weighted_frequency[[consensus_index]])
  }

  data.frame(
    alignment_col_1based = seq_len(n_columns),
    background_consensus_AA = consensus_aa,
    background_consensus_frequency = consensus_frequency,
    background_species_available = n_background,
    background_species_valid = species_valid,
    background_species_coverage = species_valid / n_background,
    background_gap_fraction = 1 - rowMeans(species_non_gap),
    stringsAsFactors = FALSE
  )
}

rolling_sum_centered <- function(values, flank = 3L) {
  n <- length(values)
  if (n == 0L) return(numeric())
  cumulative <- c(0, cumsum(values))
  centers <- seq_len(n)
  starts <- pmax(1L, centers - flank)
  ends <- pmin(n, centers + flank)
  cumulative[ends + 1L] - cumulative[starts]
}

bad_window_runs <- function(bad_window) {
  if (length(bad_window) == 0L || !any(bad_window)) {
    return(data.frame(start = integer(), end = integer(), run_length = integer()))
  }
  runs <- rle(bad_window)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  keep <- runs$values
  data.frame(start = starts[keep], end = ends[keep], run_length = runs$lengths[keep])
}

classify_global_qc <- function(match_ratio, coverage, evaluable_n) {
  if (evaluable_n < GLOBAL_MIN_EVALUABLE_CONSERVED_COLUMNS) return("insufficient_evidence")
  if (is.na(coverage) || coverage < GLOBAL_MIN_COVERAGE_WARNING) return("fail_low_coverage")
  if (is.na(match_ratio) || match_ratio < GLOBAL_MIN_MATCH_WARNING) return("fail_low_match")
  if (coverage < GLOBAL_MIN_COVERAGE_PASS) return("low_coverage_warning")
  if (match_ratio < GLOBAL_MIN_MATCH_PASS) return("low_match_warning")
  "pass"
}

collapse_reasons <- function(reasons) {
  reasons <- unique(reasons[nzchar(reasons)])
  if (length(reasons)) paste(reasons, collapse = ";") else ""
}

rare_records <- list()
trim_site_records <- list()
mafft_metric_records <- list()
protein_qc_records <- list()
eligibility_records <- list()
r_i <- t_i <- m_i <- p_i <- 0L

for (hog in hogs) {
  message("Evaluating trimAl centers and MAFFT context in ", hog)
  trim_file <- file.path(trim_dir, paste0(hog, ".trim.fa"))
  mafft_file <- file.path(mafft_dir, paste0(hog, ".mafft.fa"))
  fasta_file <- file.path(final_fasta_dir, paste0(hog, ".fa"))
  colmap_file <- file.path(colmap_dir, paste0(hog, ".colmap.csv"))

  trimmed <- readAAStringSet(trim_file)
  mafft <- readAAStringSet(mafft_file)
  original <- readAAStringSet(fasta_file)
  trim_matrix <- as.matrix(trimmed)
  mafft_matrix <- as.matrix(mafft)
  colmap <- read.csv(colmap_file, stringsAsFactors = FALSE)

  sequence_ids <- rownames(mafft_matrix)
  sequence_species <- sub("\\|.*$", "", sequence_ids)
  background_species <- setdiff(unique(sequence_species), c(TARGET_SPECIES, AUXILIARY_SPECIES))
  target_ids <- sequence_ids[sequence_species == TARGET_SPECIES]
  auxiliary_ids <- sequence_ids[sequence_species == AUXILIARY_SPECIES]
  eligible <- length(target_ids) > 0L && length(background_species) >= MIN_BACKGROUND_SPECIES

  eligibility_records[[length(eligibility_records) + 1L]] <- data.frame(
    HOG = hog,
    n_target_sequences = length(target_ids),
    n_auxiliary_sequences = length(auxiliary_ids),
    n_background_species = length(background_species),
    eligible_for_site_analysis = eligible,
    exclusion_reason = if (eligible) "" else if (length(target_ids) == 0L) "no_target_after_QC" else "fewer_than_min_background_species",
    stringsAsFactors = FALSE
  )
  if (!eligible) {
    warning("Skipping ", hog, ": insufficient target or background species.")
    next
  }

  # Layer 1: strict center-site metrics in the trimAl alignment.
  trim_species <- sub("\\|.*$", "", rownames(trim_matrix))
  trim_metrics <- compute_background_metrics(trim_matrix, trim_species, background_species)
  trim_metrics$HOG <- hog
  trim_metrics$trim_col_1based <- trim_metrics$alignment_col_1based
  trim_metrics$mafft_col_0based <- colmap$mafft_col_0based
  trim_metrics$mafft_col_1based <- colmap$mafft_col_1based
  trim_metrics$background_conserved_site <- with(
    trim_metrics,
    background_species_valid >= MIN_BACKGROUND_SPECIES &
      background_species_coverage >= MIN_BACKGROUND_SITE_COVERAGE &
      background_consensus_frequency >= MIN_BACKGROUND_CONSENSUS_FREQUENCY
  )
  trim_metrics$background_conserved_site[is.na(trim_metrics$background_conserved_site)] <- FALSE
  t_i <- t_i + 1L
  trim_site_records[[t_i]] <- trim_metrics
  write.csv(trim_metrics, file.path(site_dir, paste0(hog, ".trim_center_metrics.csv")), row.names = FALSE, na = "")

  # Layer 2: continuous context and protein-level QC in the untrimmed MAFFT alignment.
  mafft_metrics <- compute_background_metrics(mafft_matrix, sequence_species, background_species)
  mafft_metrics$HOG <- hog
  mafft_metrics$mafft_col_1based <- mafft_metrics$alignment_col_1based
  mafft_metrics$context_usable_column <- with(
    mafft_metrics,
    background_species_coverage >= CONTEXT_MIN_BACKGROUND_COVERAGE &
      background_consensus_frequency >= CONTEXT_MIN_BACKGROUND_CONSENSUS_FREQUENCY
  )
  mafft_metrics$background_conserved_column <- with(
    mafft_metrics,
    background_species_coverage >= MIN_BACKGROUND_SITE_COVERAGE &
      background_consensus_frequency >= MIN_BACKGROUND_CONSENSUS_FREQUENCY
  )
  mafft_metrics$context_usable_column[is.na(mafft_metrics$context_usable_column)] <- FALSE
  mafft_metrics$background_conserved_column[is.na(mafft_metrics$background_conserved_column)] <- FALSE
  m_i <- m_i + 1L
  mafft_metric_records[[m_i]] <- mafft_metrics

  conserved_columns <- which(mafft_metrics$background_conserved_column)
  context_usable <- mafft_metrics$context_usable_column
  context_consensus <- mafft_metrics$background_consensus_AA
  protein_states <- list()

  for (sequence_id in target_ids) {
    row <- mafft_matrix[sequence_id, ]
    evaluable_conserved <- conserved_columns[is_valid_aa(row[conserved_columns])]
    evaluable_n <- length(evaluable_conserved)
    global_match <- if (evaluable_n) mean(row[evaluable_conserved] == context_consensus[evaluable_conserved]) else NA_real_
    global_coverage <- if (length(conserved_columns)) evaluable_n / length(conserved_columns) else NA_real_
    global_qc <- classify_global_qc(global_match, global_coverage, evaluable_n)

    # Slide along original protein residue coordinates, not alignment columns.
    residue_columns <- which(row != "-")
    residue_aa <- row[residue_columns]
    residue_evaluable <- context_usable[residue_columns] & is_valid_aa(residue_aa)
    residue_matches <- residue_evaluable & residue_aa == context_consensus[residue_columns]
    window_evaluable_n <- rolling_sum_centered(as.integer(residue_evaluable), SLIDING_FLANK)
    window_match_n <- rolling_sum_centered(as.integer(residue_matches), SLIDING_FLANK)
    window_match_ratio <- ifelse(window_evaluable_n > 0, window_match_n / window_evaluable_n, NA_real_)
    evaluable_window <- window_evaluable_n >= SLIDING_MIN_EVALUABLE_COLUMNS
    bad_window <- evaluable_window & window_match_ratio < SLIDING_BAD_MATCH_THRESHOLD
    bad_n <- sum(bad_window)
    evaluable_window_n <- sum(evaluable_window)
    bad_fraction <- if (evaluable_window_n) bad_n / evaluable_window_n else NA_real_
    runs <- bad_window_runs(bad_window)
    longest_run <- if (nrow(runs)) max(runs$run_length) else 0L
    window_qc <- if (!evaluable_window_n) {
      "insufficient_evidence"
    } else if (!is.na(bad_fraction) && bad_fraction > SLIDING_FAIL_FRACTION) {
      "fail_high_bad_window_fraction"
    } else if (longest_run >= SLIDING_FAIL_RUN) {
      "fail_long_bad_window_run"
    } else if (longest_run >= SLIDING_WARNING_RUN) {
      "local_bad_window_warning"
    } else {
      "pass"
    }

    gid <- sub("^[^|]+\\|", "", sequence_id)
    p_i <- p_i + 1L
    protein_qc_records[[p_i]] <- data.frame(
      HOG = hog,
      sequence_id = sequence_id,
      GID = gid,
      protein_global_match_ratio = global_match,
      protein_global_coverage = global_coverage,
      protein_evaluable_conserved_n = evaluable_n,
      protein_background_conserved_n = length(conserved_columns),
      protein_global_qc = global_qc,
      protein_evaluable_window_n = evaluable_window_n,
      protein_bad_window_n = bad_n,
      protein_bad_window_fraction = bad_fraction,
      protein_longest_bad_window_run = longest_run,
      protein_window_qc = window_qc,
      stringsAsFactors = FALSE
    )
    protein_states[[sequence_id]] <- list(
      global_match = global_match,
      global_coverage = global_coverage,
      evaluable_conserved_n = evaluable_n,
      conserved_n = length(conserved_columns),
      global_qc = global_qc,
      evaluable_window_n = evaluable_window_n,
      bad_n = bad_n,
      bad_fraction = bad_fraction,
      longest_run = longest_run,
      window_qc = window_qc,
      bad_runs = runs
    )
  }

  strict_centers <- which(trim_metrics$background_conserved_site)
  for (trim_col in strict_centers) {
    consensus_aa <- trim_metrics$background_consensus_AA[[trim_col]]
    consensus_frequency <- trim_metrics$background_consensus_frequency[[trim_col]]
    n_valid_species <- trim_metrics$background_species_valid[[trim_col]]
    site_coverage <- trim_metrics$background_species_coverage[[trim_col]]
    mafft_col <- trim_metrics$mafft_col_1based[[trim_col]]

    target_residues <- trim_matrix[target_ids, trim_col]
    target_residues <- target_residues[is_valid_aa(target_residues)]
    target_variant_aas <- setdiff(unique(target_residues), consensus_aa)
    if (length(target_variant_aas) == 0L) next

    for (target_aa in target_variant_aas) {
      background_with_target <- 0L
      for (species in background_species) {
        residues <- trim_matrix[trim_species == species, trim_col]
        residues <- residues[is_valid_aa(residues)]
        if (length(residues) && target_aa %in% residues) background_with_target <- background_with_target + 1L
      }
      target_prevalence <- background_with_target / n_valid_species
      if (target_prevalence > MAX_BACKGROUND_TARGET_PREVALENCE) next

      smoothed_prevalence <- (background_with_target + 0.5) / (n_valid_species + 1)
      rare_score <- -log10(smoothed_prevalence)
      target_valid <- trim_matrix[target_ids, trim_col]
      target_aa_frequency <- mean(target_valid[is_valid_aa(target_valid)] == target_aa)
      auxiliary_valid <- if (length(auxiliary_ids)) trim_matrix[auxiliary_ids, trim_col] else character()
      auxiliary_valid <- auxiliary_valid[is_valid_aa(auxiliary_valid)]
      auxiliary_aa_frequency <- if (length(auxiliary_valid)) mean(auxiliary_valid == target_aa) else NA_real_
      shared_with_auxiliary <- !is.na(auxiliary_aa_frequency) && auxiliary_aa_frequency > 0
      variant_class <- if (target_aa_frequency < 1) {
        "Sbi_copy_specific"
      } else if (shared_with_auxiliary) {
        "Sbi_Sgl_shared"
      } else {
        "Sbi_fixed_not_seen_in_Sgl"
      }

      hit_ids <- target_ids[trim_matrix[target_ids, trim_col] == target_aa]
      for (sequence_id in hit_ids) {
        row <- mafft_matrix[sequence_id, ]
        gid <- sub("^[^|]+\\|", "", sequence_id)
        window_columns <- seq.int(
          max(1L, mafft_col - CONTEXT_FLANK_AA),
          min(ncol(mafft_matrix), mafft_col + CONTEXT_FLANK_AA)
        )
        flank_columns <- setdiff(window_columns, mafft_col)
        usable_flanks <- flank_columns[context_usable[flank_columns]]
        target_evaluable_flanks <- usable_flanks[is_valid_aa(row[usable_flanks])]
        target_flank_match <- if (length(target_evaluable_flanks)) {
          mean(row[target_evaluable_flanks] == context_consensus[target_evaluable_flanks])
        } else {
          NA_real_
        }

        target_gap_n <- sum(row[window_columns] == "-")
        target_gap_fraction <- target_gap_n / length(window_columns)
        background_gap_fraction <- mean(mafft_metrics$background_gap_fraction[window_columns])
        window_gap_qc <- if (length(window_columns) < MIN_WINDOW_EFFECTIVE_LENGTH) {
          "boundary_insufficient"
        } else if (target_gap_fraction > TARGET_WINDOW_MAX_GAP_FRACTION) {
          "fail_target_gap"
        } else if (background_gap_fraction > BACKGROUND_WINDOW_GAP_FAIL) {
          "fail_background_gap"
        } else if (background_gap_fraction > BACKGROUND_WINDOW_GAP_WARNING) {
          "background_gap_warning"
        } else if (target_gap_n > 0L) {
          "target_gap_allowed"
        } else {
          "pass"
        }

        protein_state <- protein_states[[sequence_id]]
        protein_position <- sum(row[seq_len(mafft_col)] != "-")
        candidate_run <- 0L
        if (nrow(protein_state$bad_runs)) {
          overlapping <- protein_state$bad_runs[
            protein_position >= pmax(1L, protein_state$bad_runs$start - SLIDING_FLANK) &
              protein_position <= protein_state$bad_runs$end + SLIDING_FLANK,
            , drop = FALSE
          ]
          if (nrow(overlapping)) candidate_run <- max(overlapping$run_length)
        }

        original_sequence <- as.character(original[[sequence_id]])
        verified <- protein_position >= 1L &&
          protein_position <= nchar(original_sequence) &&
          substr(original_sequence, protein_position, protein_position) == target_aa

        fail_reasons <- character()
        warning_reasons <- character()
        if (target_gap_fraction > TARGET_WINDOW_MAX_GAP_FRACTION) fail_reasons <- c(fail_reasons, "target_window_gap_high")
        if (background_gap_fraction > BACKGROUND_WINDOW_GAP_FAIL) {
          fail_reasons <- c(fail_reasons, "background_window_gap_high")
        } else if (background_gap_fraction > BACKGROUND_WINDOW_GAP_WARNING) {
          warning_reasons <- c(warning_reasons, "background_window_gap_warning")
        }
        if (length(window_columns) < MIN_WINDOW_EFFECTIVE_LENGTH) warning_reasons <- c(warning_reasons, "window_boundary_insufficient")
        if (length(usable_flanks) < CONTEXT_MIN_USABLE_FLANK_COLUMNS) {
          warning_reasons <- c(warning_reasons, "too_few_usable_flanks")
        } else if (length(target_evaluable_flanks) < CONTEXT_MIN_USABLE_FLANK_COLUMNS) {
          warning_reasons <- c(warning_reasons, "too_few_target_evaluable_flanks")
        } else if (is.na(target_flank_match) || target_flank_match < CONTEXT_MIN_TARGET_FLANK_MATCH) {
          fail_reasons <- c(fail_reasons, "target_flank_match_low")
        }
        if (candidate_run >= SLIDING_WARNING_RUN) fail_reasons <- c(fail_reasons, "candidate_in_consecutive_bad_windows")
        if (protein_state$global_qc %in% c("fail_low_match", "fail_low_coverage")) {
          fail_reasons <- c(fail_reasons, protein_state$global_qc)
        } else if (protein_state$global_qc != "pass") {
          warning_reasons <- c(warning_reasons, protein_state$global_qc)
        }
        if (!is.na(protein_state$bad_fraction) && protein_state$bad_fraction > SLIDING_FAIL_FRACTION) {
          fail_reasons <- c(fail_reasons, "protein_bad_window_fraction_high")
        } else if (
          protein_state$window_qc != "pass" &&
          protein_state$window_qc != "local_bad_window_warning"
        ) {
          warning_reasons <- c(warning_reasons, protein_state$window_qc)
        }
        if (!verified) fail_reasons <- c(fail_reasons, "original_coordinate_mismatch")

        confidence_class <- if (length(fail_reasons)) {
          "excluded_alignment_artifact"
        } else if (length(warning_reasons)) {
          "low_confidence"
        } else {
          "high_confidence"
        }

        r_i <- r_i + 1L
        rare_records[[r_i]] <- data.frame(
          HOG = hog, GID = gid, species = TARGET_SPECIES,
          trim_col_1based = trim_col, mafft_col_0based = mafft_col - 1L,
          mafft_col_1based = mafft_col, protein_pos_1based = protein_position,
          AA = target_aa, background_consensus_AA = consensus_aa,
          background_consensus_frequency = consensus_frequency,
          background_species_available = length(background_species),
          background_species_valid = n_valid_species,
          background_species_coverage = site_coverage,
          background_species_with_target_AA = background_with_target,
          background_target_AA_prevalence = target_prevalence,
          smoothed_background_prevalence = smoothed_prevalence,
          rare_score = rare_score, Sbi_copy_frequency = target_aa_frequency,
          Sgl_copy_frequency = auxiliary_aa_frequency,
          shared_with_Sgl = shared_with_auxiliary, variant_class = variant_class,
          window_start_mafft_1based = min(window_columns),
          window_end_mafft_1based = max(window_columns),
          window_effective_length = length(window_columns),
          flank_usable_n = length(usable_flanks),
          target_flank_evaluable_n = length(target_evaluable_flanks),
          target_flank_match_ratio = target_flank_match,
          target_window_gap_n = target_gap_n,
          target_window_gap_fraction = target_gap_fraction,
          background_window_gap_fraction = background_gap_fraction,
          window_gap_qc = window_gap_qc,
          protein_global_match_ratio = protein_state$global_match,
          protein_global_coverage = protein_state$global_coverage,
          protein_evaluable_conserved_n = protein_state$evaluable_conserved_n,
          protein_background_conserved_n = protein_state$conserved_n,
          protein_global_qc = protein_state$global_qc,
          protein_evaluable_window_n = protein_state$evaluable_window_n,
          protein_bad_window_n = protein_state$bad_n,
          protein_bad_window_fraction = protein_state$bad_fraction,
          protein_longest_bad_window_run = protein_state$longest_run,
          protein_window_qc = protein_state$window_qc,
          candidate_near_bad_window_run = candidate_run,
          original_coordinate_verified = verified,
          confidence_class = confidence_class,
          qc_fail_reasons = collapse_reasons(fail_reasons),
          qc_warning_reasons = collapse_reasons(warning_reasons),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

rare_df <- bind_rows(rare_records)
trim_site_df <- bind_rows(trim_site_records)
mafft_metric_df <- bind_rows(mafft_metric_records)
protein_qc_df <- bind_rows(protein_qc_records)
eligibility_df <- bind_rows(eligibility_records)
if (nrow(rare_df) == 0L) stop("No strict center-site candidates were detected.")

rare_hog_summary <- rare_df %>%
  group_by(HOG) %>%
  summarise(
    n_genes_with_center_candidate = n_distinct(GID),
    n_center_gene_site_records = n(),
    n_high_confidence = sum(confidence_class == "high_confidence"),
    n_low_confidence = sum(confidence_class == "low_confidence"),
    n_excluded_alignment_artifact = sum(confidence_class == "excluded_alignment_artifact"),
    n_unique_high_confidence_sites = n_distinct(
      paste(mafft_col_1based[confidence_class == "high_confidence"], AA[confidence_class == "high_confidence"])
    ),
    all_coordinates_verified = all(original_coordinate_verified),
    .groups = "drop"
  )

write.csv(rare_df, file.path(V2_TEMP, "rare_aa_records_species_weighted.csv"), row.names = FALSE, na = "")
write.csv(rare_df[rare_df$confidence_class == "high_confidence", ], file.path(V2_TEMP, "rare_aa_high_confidence.csv"), row.names = FALSE, na = "")
write.csv(rare_df[rare_df$confidence_class == "low_confidence", ], file.path(V2_TEMP, "rare_aa_low_confidence.csv"), row.names = FALSE, na = "")
write.csv(rare_df[rare_df$confidence_class == "excluded_alignment_artifact", ], file.path(V2_TEMP, "rare_aa_excluded_alignment_artifact.csv"), row.names = FALSE, na = "")
write.csv(trim_site_df, file.path(V2_TEMP, "background_site_metrics.csv"), row.names = FALSE, na = "")
write.csv(mafft_metric_df, file.path(V2_TEMP, "mafft_background_column_metrics.csv"), row.names = FALSE, na = "")
write.csv(protein_qc_df, file.path(V2_TEMP, "Sbi_protein_alignment_qc.csv"), row.names = FALSE, na = "")
write.csv(eligibility_df, file.path(V2_TEMP, "hog_analysis_eligibility.csv"), row.names = FALSE, na = "")
write.csv(rare_hog_summary, file.path(V2_TEMP, "rare_aa_hog_summary.csv"), row.names = FALSE, na = "")

saveRDS(
  list(
    parameters = list(
      center_min_background_species = MIN_BACKGROUND_SPECIES,
      center_min_background_coverage = MIN_BACKGROUND_SITE_COVERAGE,
      center_min_consensus_frequency = MIN_BACKGROUND_CONSENSUS_FREQUENCY,
      center_max_target_prevalence = MAX_BACKGROUND_TARGET_PREVALENCE,
      context_flank_aa = CONTEXT_FLANK_AA,
      context_min_background_coverage = CONTEXT_MIN_BACKGROUND_COVERAGE,
      context_min_consensus_frequency = CONTEXT_MIN_BACKGROUND_CONSENSUS_FREQUENCY,
      context_min_usable_flanks = CONTEXT_MIN_USABLE_FLANK_COLUMNS,
      context_min_target_flank_match = CONTEXT_MIN_TARGET_FLANK_MATCH,
      global_min_match_pass = GLOBAL_MIN_MATCH_PASS,
      global_min_coverage_pass = GLOBAL_MIN_COVERAGE_PASS,
      global_min_evaluable_columns = GLOBAL_MIN_EVALUABLE_CONSERVED_COLUMNS,
      sliding_window_aa = SLIDING_WINDOW_AA,
      sliding_bad_match_threshold = SLIDING_BAD_MATCH_THRESHOLD,
      target_window_max_gap_fraction = TARGET_WINDOW_MAX_GAP_FRACTION,
      background_window_gap_warning = BACKGROUND_WINDOW_GAP_WARNING,
      background_window_gap_fail = BACKGROUND_WINDOW_GAP_FAIL
    ),
    rare_aa_all_centers = rare_df,
    rare_aa_high_confidence = rare_df[rare_df$confidence_class == "high_confidence", ],
    protein_qc = protein_qc_df,
    trim_center_metrics = trim_site_df,
    mafft_background_metrics = mafft_metric_df
  ),
  file.path(V2_TEMP, "HOG_results_species_weighted.rds")
)

message(
  "Center candidates: ", nrow(rare_df),
  "; high confidence: ", sum(rare_df$confidence_class == "high_confidence"),
  "; low confidence: ", sum(rare_df$confidence_class == "low_confidence"),
  "; excluded alignment artifact: ", sum(rare_df$confidence_class == "excluded_alignment_artifact"),
  "; coordinate verification failures: ", sum(!rare_df$original_coordinate_verified)
)
