#!/usr/bin/env Rscript

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_file), "v2_pipeline_config.R"))

project_library <- file.path(PROJECT_ROOT, ".Rlib")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))
if (!requireNamespace("ggmsa", quietly = TRUE)) {
  stop(
    "The ggmsa package is required. Install it with: ",
    "dir.create('.Rlib', showWarnings=FALSE); ",
    ".libPaths(c(normalizePath('.Rlib'), .libPaths())); ",
    "BiocManager::install('ggmsa', ask=FALSE, update=FALSE)"
  )
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(ggplot2)
  library(ggmsa)
})

rare_file <- file.path(V2_TEMP, "rare_aa_records_species_weighted.csv")
mafft_dir <- file.path(V2_TEMP, "05_mafft_final")
root_candidate_file <- file.path(V2_OUTPUT, "root_candidates_v2.csv")
shoot_candidate_file <- file.path(V2_OUTPUT, "shoot_candidates_v2.csv")

scope <- Sys.getenv("RARE_AA_PLOT_SCOPE", unset = "expression_candidates")
valid_scopes <- c("expression_candidates", "all_high_confidence_sites", "all_center_candidates")
if (!scope %in% valid_scopes) {
  stop("RARE_AA_PLOT_SCOPE must be one of: ", paste(valid_scopes, collapse = ", "))
}

window_flank <- 3L
scope_directory <- switch(
  scope,
  expression_candidates = "expression_candidates",
  all_high_confidence_sites = "all_high_confidence_sites",
  all_center_candidates = "all_center_candidates"
)
plot_root <- file.path(V2_OUTPUT, "rare_AA_site_plots", scope_directory)
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

rare <- read.csv(rare_file, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(rare) == 0L) stop("No rare amino-acid records were found in ", rare_file)
if (!"confidence_class" %in% names(rare)) {
  stop("Rare-AA records lack confidence_class. Re-run 03_detect_rare_AA_v2.R first.")
}
root_candidates <- if (file.exists(root_candidate_file)) {
  read.csv(root_candidate_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(GID = character())
}
shoot_candidates <- if (file.exists(shoot_candidate_file)) {
  read.csv(shoot_candidate_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(GID = character())
}
root_candidate_ids <- unique(root_candidates$GID)
shoot_candidate_ids <- unique(shoot_candidates$GID)
expression_candidate_ids <- union(root_candidate_ids, shoot_candidate_ids)

site_key <- function(hog, column) paste(hog, column, sep = "::")
rare$key <- site_key(rare$HOG, rare$mafft_col_1based)
rare_high <- rare[rare$confidence_class == "high_confidence", , drop = FALSE]
rare_scored <- rare[
  rare$confidence_class %in% c("high_confidence", "low_confidence"),
  , drop = FALSE
]

if (scope == "expression_candidates") {
  if (length(expression_candidate_ids) == 0L) {
    stop("No root or shoot expression candidates are available for plotting.")
  }
  selected_keys <- unique(rare_scored$key[rare_scored$GID %in% expression_candidate_ids])
} else if (scope == "all_high_confidence_sites") {
  selected_keys <- unique(rare_high$key)
} else {
  selected_keys <- unique(rare$key)
}
if (length(selected_keys) == 0L) stop("No sites satisfy the selected plotting scope: ", scope)

# Once a site is selected, retain every rare Sbi protein at that HOG/column,
# even when only one of those proteins passed the expression filter.
plot_records <- rare[rare$key %in% selected_keys, , drop = FALSE]
plot_records$Root_candidate <- plot_records$GID %in% root_candidate_ids
plot_records$Shoot_candidate <- plot_records$GID %in% shoot_candidate_ids
plot_records$selected_expression_candidate <- plot_records$Root_candidate | plot_records$Shoot_candidate
plot_records$mutation_label <- paste0(
  plot_records$background_consensus_AA,
  plot_records$protein_pos_1based,
  plot_records$AA
)

site_index <- unique(plot_records[, c("HOG", "mafft_col_1based", "key")])
site_index <- site_index[order(site_index$HOG, site_index$mafft_col_1based), , drop = FALSE]

appendix_parts <- list()
plot_summary_parts <- list()
appendix_i <- summary_i <- 0L

for (hog in unique(site_index$HOG)) {
  hog_sites <- site_index[site_index$HOG == hog, , drop = FALSE]
  hog_dir <- file.path(plot_root, hog)
  dir.create(hog_dir, recursive = TRUE, showWarnings = FALSE)

  mafft_file <- file.path(mafft_dir, paste0(hog, ".mafft.fa"))
  if (!file.exists(mafft_file)) stop("Missing MAFFT alignment: ", mafft_file)
  alignment_original <- readAAStringSet(mafft_file)
  alignment_width <- unique(width(alignment_original))
  if (length(alignment_width) != 1L) stop("Sequences are not aligned in ", mafft_file)

  hog_appendix_parts <- list()
  for (site_i in seq_len(nrow(hog_sites))) {
    mafft_col <- hog_sites$mafft_col_1based[[site_i]]
    records <- plot_records[plot_records$key == hog_sites$key[[site_i]], , drop = FALSE]
    rare_gids <- unique(records$GID)
    high_confidence_gids <- unique(records$GID[records$confidence_class == "high_confidence"])
    low_confidence_gids <- unique(records$GID[records$confidence_class == "low_confidence"])
    excluded_gids <- unique(records$GID[records$confidence_class == "excluded_alignment_artifact"])

    alignment <- alignment_original
    sequence_species <- sub("\\|.*$", "", names(alignment))
    sequence_gid <- sub("^[^|]+\\|", "", names(alignment))
    sequence_confidence <- ifelse(
      sequence_gid %in% high_confidence_gids, "HC",
      ifelse(sequence_gid %in% low_confidence_gids, "LC",
             ifelse(sequence_gid %in% excluded_gids, "EXCLUDED", ""))
    )
    carries_rare_aa <- nzchar(sequence_confidence)
    species_priority <- ifelse(
      sequence_confidence == "HC", 0L,
      ifelse(sequence_confidence == "LC", 1L,
             ifelse(sequence_confidence == "EXCLUDED", 2L,
                    ifelse(sequence_species == TARGET_SPECIES, 3L,
                           ifelse(sequence_species == AUXILIARY_SPECIES, 4L, 5L))))
    )
    sequence_order <- order(species_priority, sequence_species, sequence_gid)
    alignment <- alignment[sequence_order]
    sequence_species <- sequence_species[sequence_order]
    sequence_gid <- sequence_gid[sequence_order]
    sequence_confidence <- sequence_confidence[sequence_order]
    carries_rare_aa <- carries_rare_aa[sequence_order]

    display_names <- paste(sequence_species, sequence_gid, sep = "_")
    display_names[sequence_confidence == "HC"] <- paste0("HC_RARE_", display_names[sequence_confidence == "HC"])
    display_names[sequence_confidence == "LC"] <- paste0("LC_RARE_", display_names[sequence_confidence == "LC"])
    display_names[sequence_confidence == "EXCLUDED"] <- paste0("EXCLUDED_", display_names[sequence_confidence == "EXCLUDED"])
    names(alignment) <- display_names

    plot_start <- max(1L, mafft_col - window_flank)
    plot_end <- min(alignment_width, mafft_col + window_flank)
    title_records <- records[
      records$confidence_class %in% c("high_confidence", "low_confidence"),
      , drop = FALSE
    ]
    if (nrow(title_records) == 0L) title_records <- records
    target_changes <- paste0(
      unique(title_records$background_consensus_AA),
      "->",
      paste(sort(unique(title_records$AA)), collapse = "/")
    )
    title <- paste0(hog, " | MAFFT column ", mafft_col, " | ", target_changes)
    subtitle <- paste0(
      "Window: ", plot_start, "-", plot_end,
      "; HC/LC/excluded proteins: ", length(high_confidence_gids), "/",
      length(low_confidence_gids), "/", length(excluded_gids), "\n",
      "Outlined column = candidate site; prefixes show per-protein confidence class"
    )

    p <- ggmsa(
      alignment,
      start = plot_start,
      end = plot_end,
      font = "DroidSansMono",
      color = "Chemistry_AA",
      char_width = 0.82,
      seq_name = TRUE,
      position_highlight = mafft_col,
      border = "grey55",
      show.legend = FALSE
    ) +
      labs(title = title, subtitle = subtitle) +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8),
        plot.margin = margin(8, 12, 8, 8)
      )

    plot_name <- sprintf("site_mafft_col_%06d.pdf", mafft_col)
    plot_file <- file.path(hog_dir, plot_name)
    plot_height <- max(4.5, 0.24 * length(alignment) + 2.8)
    ggsave(
      plot_file,
      p,
      width = 11,
      height = plot_height,
      units = "in",
      limitsize = FALSE
    )

    records$plot_scope <- scope
    records$plot_id <- sub("\\.pdf$", "", plot_name)
    records$plot_file <- normalizePath(plot_file, mustWork = TRUE)
    records$window_start_mafft_1based <- plot_start
    records$window_end_mafft_1based <- plot_end
    records$plot_sequence_label <- paste0(
      ifelse(records$confidence_class == "high_confidence", "HC_RARE_",
             ifelse(records$confidence_class == "low_confidence", "LC_RARE_", "EXCLUDED_")),
      records$species, "_", records$GID
    )

    appendix_i <- appendix_i + 1L
    appendix_parts[[appendix_i]] <- records
    hog_appendix_parts[[site_i]] <- records

    summary_i <- summary_i + 1L
    plot_summary_parts[[summary_i]] <- data.frame(
      plot_scope = scope,
      HOG = hog,
      plot_id = sub("\\.pdf$", "", plot_name),
      plot_file = normalizePath(plot_file, mustWork = TRUE),
      mafft_col_1based = mafft_col,
      window_start_mafft_1based = plot_start,
      window_end_mafft_1based = plot_end,
      background_consensus_AA = paste(unique(records$background_consensus_AA), collapse = "/"),
      target_AA = paste(sort(unique(records$AA)), collapse = "/"),
      n_rare_proteins = length(rare_gids),
      n_high_confidence_proteins = length(high_confidence_gids),
      n_low_confidence_proteins = length(low_confidence_gids),
      n_excluded_proteins = length(excluded_gids),
      n_root_candidate_proteins = length(intersect(rare_gids, root_candidate_ids)),
      n_shoot_candidate_proteins = length(intersect(rare_gids, shoot_candidate_ids)),
      stringsAsFactors = FALSE
    )

    if (site_i %% 25L == 0L || site_i == nrow(hog_sites)) {
      message("Plotted ", site_i, " / ", nrow(hog_sites), " sites for ", hog)
    }
  }

  hog_appendix <- bind_rows(hog_appendix_parts)
  write.csv(
    hog_appendix,
    file.path(hog_dir, paste0(hog, "_site_plot_annotations.csv")),
    row.names = FALSE,
    na = ""
  )
}

global_appendix <- bind_rows(appendix_parts) %>%
  select(
    plot_scope, HOG, plot_id, plot_file,
    window_start_mafft_1based, window_end_mafft_1based,
    GID, species, plot_sequence_label, confidence_class,
    mutation_label, qc_fail_reasons, qc_warning_reasons,
    Root_candidate, Shoot_candidate, selected_expression_candidate,
    everything()
  )
plot_summary <- bind_rows(plot_summary_parts)

write.csv(global_appendix, file.path(plot_root, "all_site_plot_annotations.csv"), row.names = FALSE, na = "")
write.csv(plot_summary, file.path(plot_root, "site_plot_summary.csv"), row.names = FALSE, na = "")

message(
  "Created ", nrow(plot_summary), " ggmsa site plots in ",
  length(unique(plot_summary$HOG)), " HOG folders. Output: ", plot_root
)
