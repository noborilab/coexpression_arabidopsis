#!/usr/bin/env Rscript
# run_pipeline.R — Config-driven CLI entrypoint for the Extended Gene Coexpression
# Analysis Pipeline.
#
# Usage:
#   Rscript inst/scripts/run_pipeline.R --config path/to/config.yaml
#   Rscript inst/scripts/run_pipeline.R --config path/to/config.yaml --validate-only

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-c", "--config"),
    type    = "character",
    default = NULL,
    help    = "Path to YAML configuration file [required]",
    metavar = "FILE"
  ),
  make_option(
    c("--validate-only"),
    action  = "store_true",
    default = FALSE,
    help    = "Validate and echo the config, then exit without running the pipeline"
  )
)

parser <- OptionParser(
  usage       = "%prog --config <config.yaml> [--validate-only]",
  option_list = option_list,
  description = paste(
    "Extended Gene Coexpression Analysis Pipeline",
    "Takes a preprocessed 10x Seurat object and produces context-robust",
    "co-expression modules with biological interpretation.",
    sep = "\n"
  )
)

args <- parse_args(parser)

if (is.null(args$config)) {
  print_help(parser)
  stop("--config is required", call. = FALSE)
}

if (!file.exists(args$config)) {
  stop(sprintf("Config file not found: %s", args$config), call. = FALSE)
}

# ---------------------------------------------------------------------------
# Load and validate config
# ---------------------------------------------------------------------------
message("Loading config: ", args$config)
cfg <- yaml::read_yaml(args$config)

required_top <- c("input", "stratum", "estimation", "output")
missing_keys <- setdiff(required_top, names(cfg))
if (length(missing_keys) > 0) {
  stop(
    "Config is missing required top-level keys: ",
    paste(missing_keys, collapse = ", "),
    call. = FALSE
  )
}

if (is.null(cfg$input$path))       stop("config$input$path is required",       call. = FALSE)
if (is.null(cfg$input$dataset_id)) stop("config$input$dataset_id is required",  call. = FALSE)
if (is.null(cfg$estimation$mode))  stop("config$estimation$mode is required",   call. = FALSE)

valid_modes <- c("pseudobulk", "singlecellggm")
if (!cfg$estimation$mode %in% valid_modes) {
  stop(sprintf(
    "config$estimation$mode must be one of: %s  (got '%s')",
    paste(valid_modes, collapse = ", "), cfg$estimation$mode
  ), call. = FALSE)
}

robustness_status <- if (isTRUE(cfg$robustness$enabled)) "enabled" else "disabled"

message("Config validated successfully.")
message("  dataset_id  : ", cfg$input$dataset_id)
message("  input path  : ", cfg$input$path)
message("  mode        : ", cfg$estimation$mode)
message("  robustness  : ", robustness_status)
message("  output dir  : ", cfg$output$dir)

if (args[["validate-only"]]) {
  message("\n--validate-only: exiting without running the pipeline.")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Pipeline execution
# TODO (Phase 2): load CoexprArabidopsis, dispatch to the appropriate
#                 adapter + estimation mode + robustness + interpretation.
# ---------------------------------------------------------------------------
stop(paste(
  "Pipeline execution not yet implemented — scaffold stub only.",
  "Use --validate-only to test your config file."
), call. = FALSE)
