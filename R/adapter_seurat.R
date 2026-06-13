#' @title Seurat Input Adapter
#'
#' @description
#' Converts a Seurat object into the core pipeline abstraction:
#' a named list of `(counts, meta, stratum_spec, dataset_id)`.
#'
#' **This is the ONLY file in the codebase that imports or depends on Seurat.**
#' All downstream core functions operate on the abstract representation, so
#' future adapters (AnnData, raw-count re-normalization, etc.) can be added
#' without touching any core logic.
#'
#' @section Input adapter contract:
#' Every adapter must return a named list with:
#' \describe{
#'   \item{`counts`}{Genes × cells normalized count matrix (sparse or dense).}
#'   \item{`meta`}{`data.frame` of cell metadata; rows = cells, columns = covariates.}
#'   \item{`stratum_spec`}{Character vector naming the stratum/context columns in `meta`.}
#'   \item{`dataset_id`}{Short string identifier for this dataset.}
#' }
#'
#' @name adapter_seurat
NULL

#' Load and validate a Seurat object, returning the core input abstraction
#'
#' @param path Path to `.rds` file containing a Seurat object.
#' @param assay Seurat assay to extract normalized counts from. Default `"RNA"`.
#' @param dataset_id Short identifier string for this dataset (used in output filenames).
#' @param stratum_vars Character vector of metadata column names that define
#'   the strata/contexts for per-context network estimation.
#'
#' @return Named list conforming to the core pipeline input contract:
#'   `counts`, `meta`, `stratum_spec`, `dataset_id`.
#' @export
load_seurat <- function(path, assay = "RNA", dataset_id, stratum_vars) {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}
