#' @title Gene-of-Interest (GOI) Lookup Resource
#'
#' @description
#' Builds and queries a per-gene lookup table for a user-supplied list of
#' genes of interest. For each GOI the table reports: module assignment,
#' module eigengene correlation (kME), hub flag, cross-context preservation
#' score, and top co-module partners.
#'
#' The GOI list path is supplied via `config$goi$path` (one AGI ID per line).
#'
#' @name goi_lookup
NULL

#' Build a GOI lookup table
#'
#' @param goi_ids Character vector of gene IDs of interest
#'   (AGI format, e.g. `"AT1G01010"`).
#' @param modules Output of [build_modules()].
#' @param preservation Optional output of [compute_preservation()].
#' @param hub_quantile kME quantile above which a gene is flagged as a hub.
#'   Default `0.95`.
#'
#' @return `data.frame` with one row per GOI and columns: `gene_id`, `module`,
#'   `kME`, `is_hub`, `preservation_zscore` (if `preservation` supplied),
#'   and `top_comodule_partners`.
#' @export
build_goi_table <- function(goi_ids,
                             modules,
                             preservation  = NULL,
                             hub_quantile  = 0.95) {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}

#' Query the GOI lookup table for a single gene
#'
#' @param gene_id A single gene ID (AGI format, e.g. `"AT1G01010"`).
#' @param goi_table Output of [build_goi_table()].
#'
#' @return Named list of lookup results for the requested gene, or `NULL`
#'   if the gene is not found in the table.
#' @export
query_goi <- function(gene_id, goi_table) {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}
