#' @title NetworkResult I/O helpers
#'
#' @description
#' Round-trip helpers for saving and reloading NetworkResult lists to/from disk.
#' Used to reconstruct networks from the per-condition edge_table.csv files
#' written by the estimation step.
#'
#' @name io
NULL

#' Load saved edge tables back into a NetworkResult list
#'
#' @param output_dir Parent directory containing per-stratum subdirectories
#'   (e.g. `"output_per_condition"`).
#' @param strata Character vector of stratum names (subdirectory names) to load.
#'   `NULL` (default) auto-detects all subdirectories.
#' @param mode Estimation mode recorded in each NetworkResult. One of
#'   `"pseudobulk"` or `"singlecellggm"`. Default `"singlecellggm"`.
#' @return Named list of NetworkResult, one element per stratum.
#' @export
load_network_results <- function(output_dir, strata = NULL,
                                 mode = "singlecellggm") {
  if (!dir.exists(output_dir))
    stop("output_dir does not exist: ", output_dir)

  if (is.null(strata)) {
    strata <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
    strata <- strata[nchar(strata) > 0L]
    if (length(strata) == 0L)
      stop("No subdirectories found in output_dir: ", output_dir)
  }

  result <- lapply(strata, function(s) {
    stratum_dir <- file.path(output_dir, s)
    edge_path   <- file.path(stratum_dir, "edge_table.csv")
    params_path <- file.path(stratum_dir, "params.json")

    if (!file.exists(edge_path))
      stop("edge_table.csv not found for stratum '", s, "': ", edge_path)

    et <- read.csv(edge_path, stringsAsFactors = FALSE)

    gene_ids <- sort(union(et$gene_id_A, et$gene_id_B))

    params <- if (file.exists(params_path)) {
      jsonlite::read_json(params_path)
    } else {
      list()
    }

    list(
      edge_table = et,
      gene_ids   = gene_ids,
      stratum_id = s,
      mode       = mode,
      params     = params,
      timestamp  = file.mtime(edge_path)
    )
  })

  setNames(result, strata)
}

#' Save a list of NetworkResults to disk
#'
#' Writes `edge_table.csv` and `params.json` into `{output_dir}/{stratum_id}/`
#' for each NetworkResult in `network_list`.
#'
#' @param network_list Named list of NetworkResult objects.
#' @param output_dir Parent directory. Per-stratum subdirectories are created
#'   as needed.
#' @return `output_dir` invisibly.
#' @export
save_network_results <- function(network_list, output_dir) {
  if (!is.list(network_list) || is.null(names(network_list)))
    stop("network_list must be a named list of NetworkResult.")

  for (nm in names(network_list)) {
    nr          <- network_list[[nm]]
    stratum_dir <- file.path(output_dir, nm)
    dir.create(stratum_dir, showWarnings = FALSE, recursive = TRUE)

    write.csv(nr$edge_table,
              file.path(stratum_dir, "edge_table.csv"),
              row.names = FALSE)

    jsonlite::write_json(nr$params,
                         file.path(stratum_dir, "params.json"),
                         auto_unbox = TRUE)
  }

  invisible(output_dir)
}
