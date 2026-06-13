#' @title Seurat Input Adapter
#'
#' @description
#' Converts a Seurat object into the pipeline's core input abstraction
#' (InputBundle): a named list of \code{(counts, cell_meta, gene_meta,
#' stratum_spec, dataset_id)}.
#'
#' \strong{This is the ONLY file in the codebase that imports or depends on
#' Seurat / SeuratObject.} All downstream core functions operate on the
#' abstract InputBundle, so future adapters (AnnData, raw-count
#' re-normalisation, etc.) can be added without touching any core logic.
#'
#' @section InputBundle contract:
#' \describe{
#'   \item{\code{counts}}{Genes \eqn{\times} cells count matrix.
#'     Row names are AT-IDs (Araport11).}
#'   \item{\code{cell_meta}}{\code{data.frame} of cell metadata; always
#'     includes a \code{cell_id} column.}
#'   \item{\code{gene_meta}}{\code{data.frame} with columns \code{gene_id}
#'     (AT-ID) and \code{gene_symbol} (NA when not recoverable).}
#'   \item{\code{stratum_spec}}{Named list: \code{$variable} (character) and
#'     \code{$levels} (character vector of ordered stratum levels).}
#'   \item{\code{dataset_id}}{Short string identifier for the dataset.}
#' }
#'
#' @name adapter_seurat
NULL

#' Load and standardise a Seurat object into an InputBundle
#'
#' @param seurat_path  Path to .rds file containing a Seurat object.
#' @param dataset_id   Short identifier string, e.g. \code{"pathogen_multiome"}.
#' @param stratum_var  Column name in \code{Seurat@@meta.data} to use for
#'   stratification (e.g. \code{"condition"}, \code{"organ"}).
#' @param stratum_levels Optional character vector of levels to keep and
#'   process, in order. \code{NULL} = use all levels found in
#'   \code{stratum_var}, sorted alphabetically.
#' @param group_var    Column name in \code{Seurat@@meta.data} for pseudobulk
#'   aggregation (e.g. \code{"sample_id"}, \code{"cluster"}). \code{NULL} if
#'   pseudobulk not used.
#' @param assay        Which Seurat assay to extract counts from.
#'   Default \code{"RNA"}.
#' @param slot         Which slot within the assay. Default \code{"data"}
#'   (log-normalised).
#' @param min_cells    Minimum number of cells expressing a gene (count \eqn{>
#'   0}) for the gene to be retained. Default \code{10}.
#' @param symbol_map   Optional \code{data.frame} with columns
#'   \code{gene_symbol} and \code{gene_id} (AT-ID, Araport11). If \code{NULL},
#'   row names are assumed to be AT-IDs already.
#' @return An InputBundle (named list; see \code{docs/OUTPUT_SCHEMA.md}).
#' @export
load_seurat <- function(seurat_path,
                        dataset_id,
                        stratum_var,
                        stratum_levels = NULL,
                        group_var      = NULL,
                        assay          = "RNA",
                        slot           = "data",
                        min_cells      = 10,
                        symbol_map     = NULL) {

  # 1. Load Seurat object
  if (!file.exists(seurat_path)) {
    stop("File not found: ", seurat_path)
  }
  obj <- readRDS(seurat_path)
  message("Loaded Seurat object: ", nrow(obj), " genes x ", ncol(obj), " cells")

  # 2. Extract count matrix (Seurat-version-agnostic via SeuratObject API)
  available_assays <- SeuratObject::Assays(obj)
  if (!assay %in% available_assays) {
    stop("Assay '", assay, "' not found in Seurat object. ",
         "Available assays: ", paste(available_assays, collapse = ", "))
  }
  # SeuratObject v5 renamed 'slot' to 'layer'; detect at runtime so the adapter
  # works with both v4 and v5 without version-pinning.
  so_ver <- tryCatch(utils::packageVersion("SeuratObject"),
                     error = function(e) package_version("0.0.0"))
  counts <- tryCatch(
    if (so_ver >= "5.0.0") {
      SeuratObject::GetAssayData(obj, assay = assay, layer = slot)
    } else {
      SeuratObject::GetAssayData(obj, assay = assay, slot = slot)
    },
    error = function(e) {
      stop("Could not extract slot/layer '", slot, "' from assay '", assay,
           "': ", conditionMessage(e))
    }
  )
  if (is.null(counts) || nrow(counts) == 0) {
    stop("Slot '", slot, "' in assay '", assay, "' returned an empty matrix.")
  }

  # 3. Gene ID handling (FLAG-04)
  #    Canonical AT-ID regex: AT + chromosome [1-5, M, C] + G + exactly 5 digits
  at_id_re  <- "^AT[1-5MC]G[0-9]{5}$"
  rn        <- rownames(counts)
  all_at_id <- all(grepl(at_id_re, rn))

  if (all_at_id) {
    gene_id     <- rn
    gene_symbol <- rep(NA_character_, length(rn))

  } else if (!is.null(symbol_map)) {
    if (!all(c("gene_symbol", "gene_id") %in% names(symbol_map))) {
      stop("symbol_map must have columns 'gene_symbol' and 'gene_id'.")
    }
    map_idx   <- match(rn, symbol_map$gene_symbol)
    matched   <- !is.na(map_idx)
    n_dropped <- sum(!matched)
    if (n_dropped > 0) {
      warning(n_dropped, " genes dropped: no AT-ID mapping found.")
    }
    gene_symbol      <- rn[matched]
    gene_id          <- symbol_map$gene_id[map_idx[matched]]
    counts           <- counts[matched, , drop = FALSE]
    rownames(counts) <- gene_id

  } else {
    warning("Row names do not look like AT-IDs and no symbol_map provided. ",
            "Proceeding with gene symbols as gene_id. ",
            "Downstream functions will fail on non-AT-IDs.")
    gene_id     <- rn
    gene_symbol <- rn
  }

  gene_meta <- data.frame(
    gene_id     = gene_id,
    gene_symbol = gene_symbol,
    stringsAsFactors = FALSE,
    row.names   = NULL
  )

  # 4. Expression filter: retain genes detected in >= min_cells cells
  #    Use Matrix::rowSums so S4 dispatch is explicit for sparse dgCMatrix objects
  #    returned by SeuratObject::GetAssayData (base rowSums is not S4-aware in
  #    package namespaces where Matrix is loaded but not attached).
  n_before    <- nrow(counts)
  n_expressed <- Matrix::rowSums(counts != 0)
  keep_genes  <- n_expressed >= min_cells
  counts    <- counts[keep_genes, , drop = FALSE]
  gene_meta <- gene_meta[keep_genes, , drop = FALSE]
  rownames(gene_meta) <- NULL
  message(sum(keep_genes), " / ", n_before,
          " genes retained after min_cells filter.")

  if (nrow(counts) == 0) {
    stop("Zero genes retained after min_cells filter (min_cells = ",
         min_cells, ").")
  }

  # 5. Cell metadata
  cell_meta         <- as.data.frame(obj@meta.data, stringsAsFactors = FALSE)
  cell_meta$cell_id <- colnames(obj)

  if (!stratum_var %in% names(cell_meta)) {
    stop("stratum_var '", stratum_var, "' not found in cell metadata. ",
         "Available columns: ", paste(names(cell_meta), collapse = ", "))
  }
  if (!is.null(group_var) && !group_var %in% names(cell_meta)) {
    stop("group_var '", group_var, "' not found in cell metadata. ",
         "Available columns: ", paste(names(cell_meta), collapse = ", "))
  }

  if (is.null(stratum_levels)) {
    stratum_levels <- sort(unique(as.character(cell_meta[[stratum_var]])))
  }

  keep_cells <- cell_meta[[stratum_var]] %in% stratum_levels
  cell_meta  <- cell_meta[keep_cells, , drop = FALSE]
  counts     <- counts[, keep_cells, drop = FALSE]

  if (ncol(counts) == 0) {
    stop("Zero cells retained after stratum_levels filtering.")
  }
  message(ncol(counts), " cells retained across ",
          length(stratum_levels), " stratum levels.")

  # 6. Assemble and return InputBundle
  list(
    counts       = counts,
    cell_meta    = cell_meta,
    gene_meta    = gene_meta,
    stratum_spec = list(variable = stratum_var, levels = stratum_levels),
    dataset_id   = dataset_id
  )
}
