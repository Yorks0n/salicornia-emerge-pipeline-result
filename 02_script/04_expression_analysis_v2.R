#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(dirname(script_file), "v2_pipeline_config.R"))

members <- read.csv(file.path(V2_TEMP, "final_member_metadata.csv"), stringsAsFactors = FALSE)
members <- members[members$species == TARGET_SPECIES, , drop = FALSE]
rare <- read.csv(file.path(V2_TEMP, "rare_aa_records_species_weighted.csv"), stringsAsFactors = FALSE)
expression <- read.csv(EXPRESSION_FILE, check.names = FALSE, stringsAsFactors = FALSE)
annotation <- read.csv(CANDIDATE_HOG_FILE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"X" %in% names(expression)) stop("Expression table must contain the gene-ID column 'X'.")
if (anyDuplicated(expression$X)) stop("Expression table contains duplicated gene IDs.")

if (!"confidence_class" %in% names(rare)) {
  stop("Rare-AA records lack confidence_class. Re-run 03_detect_rare_AA_v2.R first.")
}

rare_all_summary <- rare %>%
  group_by(GID) %>%
  summarise(
    n_center_candidates_all = n(),
    n_center_sites_all = n_distinct(paste(trim_col_1based, AA)),
    sum_rare_aa_score_all_centers = sum(rare_score),
    n_high_confidence_records = sum(confidence_class == "high_confidence"),
    n_low_confidence_records = sum(confidence_class == "low_confidence"),
    n_excluded_alignment_artifact = sum(confidence_class == "excluded_alignment_artifact"),
    .groups = "drop"
  )

rare_scored <- rare[
  rare$confidence_class %in% c("high_confidence", "low_confidence"),
  , drop = FALSE
]
rare_scored$score_weight <- ifelse(
  rare_scored$confidence_class == "high_confidence",
  1,
  LOW_CONFIDENCE_SCORE_WEIGHT
)
rare_scored$weighted_rare_score <- rare_scored$rare_score * rare_scored$score_weight

rare_scored_summary <- if (nrow(rare_scored) > 0L) {
  rare_scored %>%
    group_by(GID) %>%
    summarise(
      n_rare_AA = n(),
      n_rare_sites = n_distinct(paste(trim_col_1based, AA)),
      sum_rare_aa_score_high = sum(rare_score[confidence_class == "high_confidence"]),
      sum_rare_aa_score_low_raw = sum(rare_score[confidence_class == "low_confidence"]),
      sum_rare_aa_score_low_weighted = LOW_CONFIDENCE_SCORE_WEIGHT * sum_rare_aa_score_low_raw,
      sum_rare_aa_score = sum(weighted_rare_score),
      max_rare_aa_score_high = max(c(0, rare_score[confidence_class == "high_confidence"])),
      max_rare_aa_score_low_raw = max(c(0, rare_score[confidence_class == "low_confidence"])),
      max_rare_aa_score = max(weighted_rare_score),
      .groups = "drop"
    )
} else {
  data.frame(
    GID = character(), n_rare_AA = integer(), n_rare_sites = integer(),
    sum_rare_aa_score_high = numeric(), sum_rare_aa_score_low_raw = numeric(),
    sum_rare_aa_score_low_weighted = numeric(), sum_rare_aa_score = numeric(),
    max_rare_aa_score_high = numeric(), max_rare_aa_score_low_raw = numeric(),
    max_rare_aa_score = numeric()
  )
}
rare_summary <- full_join(rare_all_summary, rare_scored_summary, by = "GID")

detail <- members %>%
  select(HOG, OG, GID, sequence_id, sequence_length, length_ratio, trimmed_alignment_coverage) %>%
  left_join(expression, by = c("GID" = "X")) %>%
  left_join(rare_summary, by = "GID")

rare_numeric_columns <- c(
  "n_center_candidates_all", "n_center_sites_all", "sum_rare_aa_score_all_centers",
  "n_high_confidence_records", "n_low_confidence_records",
  "n_excluded_alignment_artifact", "n_rare_AA", "n_rare_sites",
  "sum_rare_aa_score_high", "sum_rare_aa_score_low_raw",
  "sum_rare_aa_score_low_weighted", "sum_rare_aa_score",
  "max_rare_aa_score_high", "max_rare_aa_score_low_raw", "max_rare_aa_score"
)
for (column in rare_numeric_columns) detail[[column]][is.na(detail[[column]])] <- 0
detail$low_confidence_score_weight <- LOW_CONFIDENCE_SCORE_WEIGHT
detail$rare_score_basis <- sprintf(
  "high_confidence_plus_%.2f_x_low_confidence",
  LOW_CONFIDENCE_SCORE_WEIGHT
)

expression_status <- function(gid, value, expression_ids) {
  if (!(gid %in% expression_ids)) return("missing_id")
  if (is.na(value)) return("missing_measurement")
  if (value == 0) return("observed_zero")
  if (value > 0) return("observed_positive")
  "invalid_negative"
}

add_tissue_metrics <- function(data, tissue_column, prefix) {
  values <- data[[tissue_column]]
  status <- mapply(expression_status, data$GID, values, MoreArgs = list(expression_ids = expression$X))
  data[[paste0(prefix, "_status")]] <- status
  valid <- status %in% c("observed_positive", "observed_zero")

  data[[paste0(prefix, "_HOG_n_valid")]] <- ave(valid, data$HOG, FUN = sum)
  data[[paste0(prefix, "_HOG_n_missing")]] <- ave(!valid, data$HOG, FUN = sum)
  data[[paste0(prefix, "_HOG_median_TPM")]] <- ave(
    ifelse(valid, values, NA_real_), data$HOG,
    FUN = function(x) if (sum(!is.na(x)) >= MIN_VALID_EXPRESSION_MEMBERS) median(x, na.rm = TRUE) else NA_real_
  )
  data[[paste0(prefix, "_HOG_total_TPM")]] <- ave(
    ifelse(valid, values, NA_real_), data$HOG,
    FUN = function(x) if (sum(!is.na(x)) >= MIN_VALID_EXPRESSION_MEMBERS) sum(x, na.rm = TRUE) else NA_real_
  )
  data[[paste0(prefix, "_HOG_rank")]] <- ave(
    ifelse(valid, values, NA_real_), data$HOG,
    FUN = function(x) rank(-x, ties.method = "min", na.last = "keep")
  )
  total <- data[[paste0(prefix, "_HOG_total_TPM")]]
  data[[paste0(prefix, "_expression_share")]] <- ifelse(valid & total > 0, values / total, NA_real_)

  median_baseline <- data[[paste0(prefix, "_HOG_median_TPM")]]
  data[[paste0(prefix, "_log2_ratio_to_median_p0p1")]] <- ifelse(
    valid & !is.na(median_baseline),
    log2((values + EXPRESSION_PSEUDOCOUNT) / (median_baseline + EXPRESSION_PSEUDOCOUNT)),
    NA_real_
  )
  data
}

detail <- add_tissue_metrics(detail, "Root_salt", "Root")
detail <- add_tissue_metrics(detail, "Shoot_salt", "Shoot")

annotation_columns <- intersect(
  c("HOG", "SYMBOL", "SYMBOL_NAME", "DESCRIPTION.x", "DESCRIPTION.y", "Preferred_name"),
  names(annotation)
)
annotation_small <- annotation[, annotation_columns, drop = FALSE]
annotation_small <- annotation_small[!duplicated(annotation_small$HOG), , drop = FALSE]
detail <- left_join(detail, annotation_small, by = "HOG")

root_score_col <- "Root_log2_ratio_to_median_p0p1"
shoot_score_col <- "Shoot_log2_ratio_to_median_p0p1"

root_score <- detail[[root_score_col]]
shoot_score <- detail[[shoot_score_col]]
detail$Root_candidate <- detail$Root_HOG_n_valid >= MIN_VALID_EXPRESSION_MEMBERS &
  detail$sum_rare_aa_score > MIN_RARE_SCORE_SUM &
  !is.na(root_score) & root_score > MIN_EXPRESSION_LOG2_RATIO
detail$Shoot_candidate <- detail$Shoot_HOG_n_valid >= MIN_VALID_EXPRESSION_MEMBERS &
  detail$sum_rare_aa_score > MIN_RARE_SCORE_SUM &
  !is.na(shoot_score) & shoot_score > MIN_EXPRESSION_LOG2_RATIO

write.csv(detail, file.path(V2_OUTPUT, "sbi_HOG_rare_detail_v2.csv"), row.names = FALSE, na = "")
write.csv(detail[detail$Root_candidate, , drop = FALSE], file.path(V2_OUTPUT, "root_candidates_v2.csv"), row.names = FALSE, na = "")
write.csv(detail[detail$Shoot_candidate, , drop = FALSE], file.path(V2_OUTPUT, "shoot_candidates_v2.csv"), row.names = FALSE, na = "")

status_count <- function(status, tissue) {
  counts <- as.data.frame(table(status = status), stringsAsFactors = FALSE)
  names(counts)[[2]] <- "n_genes"
  counts$tissue <- tissue
  counts[, c("tissue", "status", "n_genes")]
}
missing_summary <- bind_rows(
  status_count(detail$Root_status, "Root"),
  status_count(detail$Shoot_status, "Shoot")
)
write.csv(missing_summary, file.path(V2_OUTPUT, "expression_missingness_summary.csv"), row.names = FALSE)

plot_tissue <- function(data, score_column, candidate_column, x_label, output_file) {
  plot_data <- data[is.finite(data[[score_column]]) & is.finite(log2(data$sum_rare_aa_score + 1)), , drop = FALSE]
  p <- ggplot(plot_data, aes(x = .data[[score_column]], y = log2(sum_rare_aa_score + 1))) +
    geom_point(aes(color = .data[[candidate_column]]), alpha = 0.55) +
    scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#C51B3A"), guide = "none") +
    geom_vline(xintercept = MIN_EXPRESSION_LOG2_RATIO, linetype = "dashed", color = "#C51B3A") +
    geom_hline(yintercept = log2(MIN_RARE_SCORE_SUM + 1), linetype = "dashed", color = "#C51B3A") +
    labs(x = x_label, y = "log2(weighted rare-AA score + 1)") +
    theme_classic()
  ggsave(file.path(V2_OUTPUT, output_file), p, width = 4.6, height = 4.2)
}

plot_tissue(detail, root_score_col, "Root_candidate", "Root log2((TPM + 0.1) / (HOG median TPM + 0.1))", "Root_rare_gene_v2.pdf")
plot_tissue(detail, shoot_score_col, "Shoot_candidate", "Shoot log2((TPM + 0.1) / (HOG median TPM + 0.1))", "Shoot_rare_gene_v2.pdf")

message(
  "Expression integration complete. Root candidates: ", sum(detail$Root_candidate, na.rm = TRUE),
  "; shoot candidates: ", sum(detail$Shoot_candidate, na.rm = TRUE),
  ". Outputs: ", V2_OUTPUT
)
