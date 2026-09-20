# Shared configuration for the HOG-direct rare-amino-acid pipeline.

get_top_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Run this file through Rscript so that the project path can be resolved.")
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

SCRIPT_DIR <- dirname(get_top_script_path())
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

find_program <- function(env_name, program, fallback = "") {
  configured <- Sys.getenv(env_name, unset = "")
  if (nzchar(configured) && file.exists(configured)) return(configured)
  on_path <- Sys.which(program)
  if (nzchar(on_path)) return(unname(on_path))
  if (nzchar(fallback) && file.exists(fallback)) return(fallback)
  stop("Cannot find ", program, ". Set environment variable ", env_name, ".")
}

ORTHOFINDER_ROOT <- env_or_default(
  "RARE_AA_ORTHOFINDER_ROOT",
  file.path(PROJECT_ROOT, "01_input", "orthofinder")
)
N0_FILE <- file.path(ORTHOFINDER_ROOT, "N0.tsv")
PROTEIN_DIR <- file.path(ORTHOFINDER_ROOT, "01_pep")
EXPRESSION_FILE <- env_or_default(
  "RARE_AA_EXPRESSION_FILE",
  file.path(PROJECT_ROOT, "01_input", "expression", "Sbi_TPM_mean.csv")
)
CANDIDATE_HOG_FILE <- file.path(PROJECT_ROOT, "01_input", "candidate_HOGs_with_annotation.csv")
V2_TEMP <- file.path(PROJECT_ROOT, "03_tempOut", "v2_hog_direct")
V2_OUTPUT <- file.path(PROJECT_ROOT, "04_output", "v2_hog_direct")

TARGET_SPECIES <- "Sbi"
AUXILIARY_SPECIES <- "Sgl"
METADATA_COLUMNS <- c("HOG", "OG", "Gene Tree Parent Clade")
STANDARD_AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1]]

PRE_QC_MIN_LENGTH_RATIO <- 0.5
PRE_QC_MAX_LENGTH_RATIO <- 2.0
PRE_QC_MAX_X_FRACTION <- 0.05
MIN_ALIGNMENT_COVERAGE <- 0.5

MIN_BACKGROUND_SPECIES <- 4L
MIN_BACKGROUND_SITE_COVERAGE <- 0.70
MIN_BACKGROUND_CONSENSUS_FREQUENCY <- 0.80
MAX_BACKGROUND_TARGET_PREVALENCE <- 0.05

# MAFFT-context QC. The strict center is detected in trimAl, whereas flanks,
# gaps, global similarity and sliding windows are evaluated in the untrimmed
# MAFFT alignment so that deleted columns are never concatenated artificially.
CONTEXT_FLANK_AA <- 3L
CONTEXT_MIN_BACKGROUND_COVERAGE <- 0.50
CONTEXT_MIN_BACKGROUND_CONSENSUS_FREQUENCY <- 0.60
CONTEXT_MIN_USABLE_FLANK_COLUMNS <- 3L
CONTEXT_MIN_TARGET_FLANK_MATCH <- 0.50

GLOBAL_MIN_MATCH_PASS <- 0.60
GLOBAL_MIN_MATCH_WARNING <- 0.45
GLOBAL_MIN_COVERAGE_PASS <- 0.60
GLOBAL_MIN_COVERAGE_WARNING <- 0.40
GLOBAL_MIN_EVALUABLE_CONSERVED_COLUMNS <- 30L

SLIDING_WINDOW_AA <- 7L
SLIDING_MIN_EVALUABLE_COLUMNS <- 4L
SLIDING_BAD_MATCH_THRESHOLD <- 0.40
SLIDING_WARNING_RUN <- 3L
SLIDING_FAIL_RUN <- 5L
SLIDING_FAIL_FRACTION <- 0.25

TARGET_WINDOW_MAX_GAP_FRACTION <- 0.40
BACKGROUND_WINDOW_GAP_WARNING <- 0.30
BACKGROUND_WINDOW_GAP_FAIL <- 0.45
MIN_WINDOW_EFFECTIVE_LENGTH <- 5L

EXPRESSION_PSEUDOCOUNT <- 0.1
MIN_VALID_EXPRESSION_MEMBERS <- 2L
LOW_CONFIDENCE_SCORE_WEIGHT <- 0.50
MIN_RARE_SCORE_SUM <- 10
MIN_EXPRESSION_LOG2_RATIO <- 2

MAFFT_THREADS <- as.integer(env_or_default("RARE_AA_THREADS", "4"))
MAFFT_BIN <- find_program(
  "RARE_AA_MAFFT",
  "mafft"
)
TRIMAL_BIN <- find_program(
  "RARE_AA_TRIMAL",
  "trimal"
)

dir.create(V2_TEMP, recursive = TRUE, showWarnings = FALSE)
dir.create(V2_OUTPUT, recursive = TRUE, showWarnings = FALSE)
