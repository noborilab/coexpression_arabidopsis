#' @title Pseudobulk Network Estimation
#'
#' @description
#' Aggregates cells to pseudobulk replicates per stratum, then estimates a
#' marginal Spearman correlation network. Produces one network per unique
#' combination of stratum variable values.
#'
#' Robustness statistics (R_score) are computed downstream in `R/robustness.R`.
#'
#' @name estimate_pseudobulk
NULL

#' Aggregate cells to pseudobulk replicates per stratum
#'
#' @param input Core pipeline input object returned by an adapter function
#'   (see `R/adapter_seurat.R` for the contract).
#' @param group_by Name of the metadata column used to aggregate cells into
#'   pseudobulk replicates (e.g. `"sample"`).
#' @param min_cells Minimum number of cells required per group; groups below
#'   this threshold are silently dropped. Default `10`.
#'
#' @return A named list of pseudobulk count matrices (genes × replicates),
#'   one entry per stratum.
#' @export
aggregate_pseudobulk <- function(input, group_by, min_cells = 10) {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}

#' Estimate per-stratum marginal Spearman correlation network from pseudobulk data
#'
#' @param pseudobulk Output of [aggregate_pseudobulk()].
#'
#' @return A named list of gene × gene correlation matrices, one per stratum.
#' @export
estimate_network_pseudobulk <- function(pseudobulk) {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}
