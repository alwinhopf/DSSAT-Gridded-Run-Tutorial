# =============================================================================
# run_all.R — execute the whole analysis workflow in one go.
#
# Chains the three stages so you don't have to invoke them separately:
#   1. SWEEP      run_experiment.R           — runs the weather x soil experiment
#   2. ANALYSIS   analyze_experiment.R       — maps / variance / pairwise tables
#   3. VALIDATION validate_against_observed.R (optional) — observed vs sim
#
# USAGE
#   Rscript run_all.R                         # sweep + analysis, using experiment.yml
#   Rscript run_all.R --experiment my.yml     # use a different experiment file
#   Rscript run_all.R --validate              # also run the observed-data validation
#   Rscript run_all.R --validate-only         # only the validation stage
#   Rscript run_all.R --no-analysis           # sweep only
# =============================================================================

suppressWarnings(suppressMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml", repos = "https://cloud.r-project.org")
}))

PROJECT_ROOT <- tryCatch({
  if (requireNamespace("this.path", quietly = TRUE)) dirname(this.path::this.path()) else getwd()
}, error = function(e) getwd())
setwd(normalizePath(PROJECT_ROOT, mustWork = FALSE))

args     <- commandArgs(trailingOnly = TRUE)
flag     <- function(name, default = NA) { i <- which(args == name); if (length(i) && i < length(args)) args[i + 1] else default }
has      <- function(name) any(args == name)
EXP_YML  <- flag("--experiment", "experiment.yml")
DO_SWEEP <- !has("--validate-only")
DO_ANALYSIS <- !has("--no-analysis") && !has("--validate-only")
DO_VALID <- has("--validate") || has("--validate-only")

step <- function(label, cmd, cmd_args) {
  cat(sprintf("\n========================================================\n>> %s\n   Rscript %s\n========================================================\n",
              label, paste(c(cmd, cmd_args), collapse = " ")))
  st <- system2("Rscript", args = c(shQuote(cmd), cmd_args))
  if (st != 0) stop(sprintf("Stage failed (exit %d): %s", st, label))
  invisible(st)
}

# ---- 1. Sweep -------------------------------------------------------------
if (DO_SWEEP) {
  if (!file.exists(EXP_YML)) stop("Experiment file not found: ", EXP_YML)
  step("STAGE 1/3  Weather x soil sweep", "run_experiment.R", shQuote(EXP_YML))
}

# ---- 2. Analysis ----------------------------------------------------------
if (DO_ANALYSIS) {
  cfg  <- yaml::read_yaml(EXP_YML)
  name <- if (is.null(cfg$experiment_name)) "experiment" else cfg$experiment_name
  combined <- if (!is.null(cfg$output$combined_csv)) cfg$output$combined_csv
              else file.path("results", paste0("EXPERIMENT_", name, "_combined.csv"))
  rvs <- if (!is.null(cfg$options$response_vars)) as.character(unlist(cfg$options$response_vars))
         else if (!is.null(cfg$options$response_var)) as.character(cfg$options$response_var)
         else "final_grain_kg_ha"
  if (!file.exists(combined)) {
    cat("\n[skip analysis] no combined output at ", combined, " — did the sweep produce results?\n", sep = "")
  } else for (rv in rvs) step(sprintf("STAGE 2/3  Analysis (%s)", rv),
                              "analyze_experiment.R", c(shQuote(combined), shQuote(rv)))
}

# ---- 3. Validation (optional) ---------------------------------------------
if (DO_VALID) step("STAGE 3/3  Validation vs observed data", "validate_against_observed.R", character(0))

cat("\nAll requested stages complete.\n")
