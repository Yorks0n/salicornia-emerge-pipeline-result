#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"

Rscript "$script_dir/01_prepare_HOG_members_v2.R"
Rscript "$script_dir/02_align_HOGs_v2.R"
Rscript "$script_dir/03_detect_rare_AA_v2.R"
Rscript "$script_dir/04_expression_analysis_v2.R"
Rscript "$script_dir/05_plot_rare_AA_sites_v2.R"
