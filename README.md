# Salicornia EMERGE Analysis Results

## Directory structure

```text
00_publish/
├── 01_input/
│   ├── candidate_HOGs_with_annotation.csv
│   ├── expression/Sbi_TPM_mean.csv
│   ├── orthofinder/N0.tsv
│   ├── orthofinder/01_pep/<species>.fa
│   └── references/ecoli_NarK.faa
├── 02_script/
├── 03_tempOut/
└── 04_output/
```

## Pre-run configuration

### Software and R packages

- Install R and ensure that `Rscript` is available from the command line. The workflow invokes `Rscript`, not RStudio. On macOS/Linux, check this with `command -v Rscript`. If it is not found, add R's `bin` directory to `PATH` and open a new terminal.
- Install MAFFT and trimAl. They must be executable binaries compatible with the current operating system and CPU architecture. Having only the installation directory is not sufficient: the programs must either be available through `PATH` or be specified by their executable paths below.
- Install and verify the following packages in the R environment used to run this workflow: `Biostrings`, `dplyr`, `ggplot2`, and `ggmsa`. If you switch to a different R/Rscript environment, install the packages in that environment as well.

Use the following command to check whether the current `Rscript` can find the required R packages:

```bash
Rscript -e 'pkgs <- c("Biostrings", "dplyr", "ggplot2", "ggmsa"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))'
```

### Software paths, data paths, and runtime parameters

The default configuration is defined in [`02_script/v2_pipeline_config.R`](02_script/v2_pipeline_config.R). MAFFT and trimAl paths are first read from environment variables. If these variables are not set, or if the specified files do not exist, the configuration searches for the programs by name in `PATH`. Environment variables are recommended so that personal absolute paths are not committed to the GitHub repository.

| Variable | Purpose | Default/value format |
|---|---|---|
| `RARE_AA_MAFFT` | Full path to the MAFFT executable | If unset, searches `PATH` for `mafft` |
| `RARE_AA_TRIMAL` | Full path to the trimAl executable | If unset, searches `PATH` for `trimal` |
| `RARE_AA_ORTHOFINDER_ROOT` | Directory containing both `N0.tsv` and `01_pep/` | `01_input/orthofinder/`; normally no change is needed |
| `RARE_AA_EXPRESSION_FILE` | Full path to the expression matrix | `01_input/expression/Sbi_TPM_mean.csv`; normally no change is needed |
| `RARE_AA_THREADS` | Number of threads used by MAFFT | `4` |
| `RARE_AA_PLOT_SCOPE` | Plotting scope for step 05 | `expression_candidates`; alternatives are `all_high_confidence_sites` and `all_center_candidates` |

On macOS/Linux, specify the executable paths and thread count before running the main workflow:

```bash
export RARE_AA_MAFFT="/absolute/path/to/mafft"
export RARE_AA_TRIMAL="/absolute/path/to/trimal"
export RARE_AA_THREADS=4
bash 02_script/run_v2_pipeline.sh
```

Alternatively, add the directories containing MAFFT and trimAl to `PATH` and leave the first two variables unset:

```bash
export PATH="/path/to/mafft/bin:/path/to/trimal/bin:$PATH"
command -v Rscript
command -v mafft
command -v trimal
bash 02_script/run_v2_pipeline.sh
```

The analysis thresholds and default parameters are also defined in `02_script/v2_pipeline_config.R`. If you need to change a scientific filtering threshold, first check the method description and then modify the corresponding constant in that file. Software paths should preferably be configured through the environment variables described above.

## Running the workflow

After configuring the software, R packages, and input files, run the complete workflow from the `06_publish` directory:

```bash
bash 02_script/run_v2_pipeline.sh
```

All intermediate files are written to `03_tempOut/v2_hog_direct/`, and final results are written to `04_output/v2_hog_direct/`.
