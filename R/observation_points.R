#' @title Observation Point Generators
#'
#' @description
#' Functions to aggregate single-cell count matrices into "observation points"
#' suitable for gene-gene co-expression estimation.  Each generator returns an
#' \strong{ObsPointSet}, a named list with the following slots:
#'
#' \describe{
#'   \item{\code{$matrix}}{Numeric matrix, genes × observation_points (aggregated
#'     expression profile).  Row names are AT-IDs; column names are point IDs.}
#'   \item{\code{$n_cells}}{Named integer vector, one entry per observation point,
#'     giving the number of cells that were aggregated into that point.}
#'   \item{\code{$point_meta}}{Data frame with one row per observation point.
#'     Always contains \code{point_id} (character, unique); may contain additional
#'     composition columns (dominant stratum/condition, source cluster, etc.)
#'     depending on the generator.}
#'   \item{\code{$gene_ids}}{Character vector of AT-IDs, identical to
#'     \code{rownames($matrix)}.}
#'   \item{\code{$design}}{Named list: \code{name} (character, generator name)
#'     plus all parameters used for that generator.}
#'   \item{\code{$aggregation}}{\code{"sum"} or \code{"mean"}.}
#' }
#'
#' @name ObsPointSet
NULL


# ---- .make_obs ----
# Validate raw pieces and assemble an ObsPointSet.
.make_obs <- function(matrix, n_cells, point_meta,
                      gene_ids, design_name, design_params,
                      aggregation) {
  stopifnot(
    is.matrix(matrix),
    nrow(matrix) > 0L,
    ncol(matrix) > 0L,
    identical(rownames(matrix), gene_ids),
    is.integer(n_cells),
    !is.null(names(n_cells)),
    length(n_cells) == ncol(matrix),
    is.data.frame(point_meta),
    nrow(point_meta) == ncol(matrix),
    "point_id" %in% names(point_meta),
    aggregation %in% c("sum", "mean")
  )
  list(
    matrix      = matrix,
    n_cells     = n_cells,
    point_meta  = point_meta,
    gene_ids    = gene_ids,
    design      = c(list(name = design_name), design_params),
    aggregation = aggregation
  )
}


# ---- .aggregate_cells ----
# Aggregate a genes×cells matrix by group labels; returns a dense genes×groups matrix.
.aggregate_cells <- function(mat, groups, aggregation = "mean") {
  stopifnot(aggregation %in% c("sum", "mean"))
  lvls <- unique(groups)
  out  <- vapply(lvls, function(g) {
    ci <- groups == g
    sub <- mat[, ci, drop = FALSE]
    if (aggregation == "sum") {
      as.numeric(Matrix::rowSums(sub))
    } else {
      as.numeric(Matrix::rowMeans(sub))
    }
  }, numeric(nrow(mat)))
  # vapply over multiple groups returns genes×groups; single group returns a vector
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = nrow(mat), ncol = 1L)
  }
  rownames(out) <- rownames(mat)
  colnames(out) <- lvls
  out
}


# ---- .pca_coords ----
# Compute PCA on transposed log-norm matrix; returns cells×n_pcs score matrix.
.pca_coords <- function(mat) {
  n_pcs <- min(30L, nrow(mat) - 1L, ncol(mat) - 1L)
  if (n_pcs < 1L) stop("Matrix too small for PCA (need > 1 gene and > 1 cell).")
  mat_dense <- as.matrix(mat)
  pca <- suppressOutput(
    prcomp(t(mat_dense), center = TRUE, scale. = FALSE, rank. = n_pcs)
  )
  pca$x  # cells × n_pcs
}

# Helper to suppress cat/print output from prcomp (not normally needed but safe).
suppressOutput <- function(expr) {
  capture.output(res <- expr)
  res
}


# ---- .build_knn_graph ----
# Build an undirected kNN igraph from a cells×n_pcs PCA matrix using chunked distances.
.build_knn_graph <- function(pca_mat, k) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for kNN graph construction. ",
         "Install it with: install.packages('igraph')")
  }
  n_cells    <- nrow(pca_mat)
  k_use      <- min(k, n_cells - 1L)
  norms      <- rowSums(pca_mat^2)           # n_cells precomputed squared norms
  chunk_size <- min(500L, n_cells)
  n_chunks   <- ceiling(n_cells / chunk_size)

  from_vec <- integer(0L)
  to_vec   <- integer(0L)

  for (ci in seq_len(n_chunks)) {
    idx_start <- (ci - 1L) * chunk_size + 1L
    idx_end   <- min(ci * chunk_size, n_cells)
    batch     <- pca_mat[idx_start:idx_end, , drop = FALSE]

    # squared Euclidean: ||a - b||^2 = ||a||^2 + ||b||^2 - 2 a·b
    sq_dist <- outer(rowSums(batch^2), norms, "+") - 2 * tcrossprod(batch, pca_mat)
    sq_dist <- pmax(sq_dist, 0)  # clamp numerical noise

    batch_size <- nrow(batch)
    for (li in seq_len(batch_size)) {
      global_i     <- idx_start + li - 1L
      dists_i      <- sq_dist[li, ]
      dists_i[global_i] <- Inf  # exclude self
      nn_idx       <- order(dists_i)[seq_len(k_use)]
      from_vec     <- c(from_vec, rep(global_i, k_use))
      to_vec       <- c(to_vec, nn_idx)
    }
  }

  edges_df <- data.frame(from = from_vec, to = to_vec)
  g <- igraph::graph_from_data_frame(edges_df,
                                     directed  = FALSE,
                                     vertices  = data.frame(name = seq_len(n_cells)))
  igraph::simplify(g)
}


# ---- .greedy_farthest_points ----
# Greedy farthest-point sampling; returns integer indices of n_anchors cells.
.greedy_farthest_points <- function(pca_mat, n_anchors) {
  n_cells <- nrow(pca_mat)
  n_anchors <- min(n_anchors, n_cells)

  center    <- colMeans(pca_mat)
  d_to_ctr  <- rowSums(sweep(pca_mat, 2L, center, "-")^2)
  chosen    <- which.min(d_to_ctr)  # start closest to mean

  # min distance from each cell to the chosen set
  min_dist <- rowSums(sweep(pca_mat, 2L, pca_mat[chosen, ], "-")^2)
  min_dist[chosen] <- -Inf  # never re-select

  for (i in seq_len(n_anchors - 1L)) {
    next_pt  <- which.max(min_dist)
    chosen   <- c(chosen, next_pt)
    new_dists <- rowSums(sweep(pca_mat, 2L, pca_mat[next_pt, ], "-")^2)
    min_dist  <- pmin(min_dist, new_dists)
    min_dist[chosen] <- -Inf
  }

  chosen
}


# ---------------------------------------------------------------------------
# Exported generators
# ---------------------------------------------------------------------------

#' Cluster-based observation points
#'
#' @description
#' Aggregates cells by graph-clustering labels at a given Louvain resolution.
#' The function checks \code{bundle$cell_meta} for an exact resolution-specific
#' column: \code{paste0("RNA_snn_res.", resolution)},
#' \code{paste0("SCT_snn_res.", resolution)}, or
#' \code{paste0("wsnn_res.", resolution)}.  Generic columns such as
#' \code{"seurat_clusters"} or any column whose name merely contains
#' \code{"cluster"} are \strong{never} used as a fallback, as they may reflect
#' a different resolution and make the sweep uninformative.
#' If no resolution-specific column exists the function always recomputes:
#' it builds a kNN graph (\code{k = 15}) from PCA coordinates and runs
#' \code{igraph::cluster_louvain} at the requested resolution.  The resulting
#' column is recorded as \code{RNA_snn_res.\{resolution\}} in
#' \code{$design$cluster_col}.
#'
#' @param bundle      An \code{InputBundle} as returned by \code{load_seurat()}.
#' @param resolution  Numeric; Louvain resolution controlling cluster granularity
#'   (higher → more clusters).  Default \code{1.0}.
#' @param aggregation \code{"sum"} or \code{"mean"}.  For \code{"sum"},
#'   \code{bundle$counts_raw} is used when available; otherwise a warning is
#'   emitted and \code{bundle$counts} (log-normalised) is used.
#'   For \code{"mean"}, \code{bundle$counts} is always used.  Default
#'   \code{"mean"}.
#'
#' @return An \code{ObsPointSet} (see \code{\link{ObsPointSet}}).
#' @export
obs_cluster <- function(bundle, resolution = 1.0, aggregation = "mean") {
  tryCatch({
    stopifnot(
      is.list(bundle),
      is.matrix(bundle$counts) || inherits(bundle$counts, "Matrix"),
      is.data.frame(bundle$cell_meta),
      aggregation %in% c("sum", "mean")
    )

    cell_meta <- bundle$cell_meta

    # --- select the count matrix ---
    if (aggregation == "sum") {
      if (is.null(bundle$counts_raw)) {
        warning("obs_cluster: counts_raw not available; falling back to ",
                "sum of log-normalized counts.")
        count_mat <- bundle$counts
      } else {
        count_mat <- bundle$counts_raw
      }
    } else {
      count_mat <- bundle$counts
    }

    # --- find or compute cluster labels ---
    # Only resolution-specific columns are accepted as reuse candidates.
    # Generic columns ("seurat_clusters", grep on "cluster") are intentionally
    # excluded: they may reflect a different resolution and corrupt a granularity
    # sweep (FLAG-14 Bug #1).
    cluster_col <- NULL
    candidates <- c(
      paste0("RNA_snn_res.", resolution),
      paste0("SCT_snn_res.", resolution),
      paste0("wsnn_res.", resolution)
    )
    for (cand in candidates) {
      if (cand %in% names(cell_meta)) {
        cluster_col <- cand
        break
      }
    }

    if (!is.null(cluster_col)) {
      labels <- as.character(cell_meta[[cluster_col]])
    } else {
      if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required for Louvain clustering. ",
             "Install it with: install.packages('igraph')")
      }
      pca_mat <- .pca_coords(bundle$counts)
      g       <- .build_knn_graph(pca_mat, k = 15L)
      comm    <- tryCatch(
        igraph::cluster_louvain(g, resolution = resolution),
        error = function(e) igraph::cluster_louvain(g)
      )
      labels      <- as.character(igraph::membership(comm))
      cluster_col <- paste0("RNA_snn_res.", resolution)
    }

    # --- aggregate ---
    agg_mat <- .aggregate_cells(count_mat, labels, aggregation)
    agg_mat <- as.matrix(agg_mat)

    lvls      <- colnames(agg_mat)
    n_per     <- setNames(
      as.integer(table(labels)[lvls]),
      lvls
    )
    point_ids <- paste0("cluster_", lvls)
    colnames(agg_mat) <- point_ids
    names(n_per)      <- point_ids

    point_meta <- data.frame(
      point_id       = point_ids,
      source_cluster = lvls,
      stringsAsFactors = FALSE
    )

    message("obs_cluster: ", length(point_ids), " points, median ",
            median(n_per), " cells/point")

    .make_obs(
      matrix      = agg_mat,
      n_cells     = n_per,
      point_meta  = point_meta,
      gene_ids    = rownames(agg_mat),
      design_name = "obs_cluster",
      design_params = list(
        resolution  = resolution,
        cluster_col = cluster_col
      ),
      aggregation = aggregation
    )
  }, error = function(e) {
    stop("obs_cluster failed: ", conditionMessage(e))
  })
}


#' Subcluster-based observation points
#'
#' @description
#' Aggregates cells by a precomputed subcluster column already present in
#' \code{bundle$cell_meta}.  Each unique value of \code{group_col} becomes one
#' observation point.
#'
#' @param bundle      An \code{InputBundle}.
#' @param group_col   Name of the column in \code{bundle$cell_meta} that holds
#'   the subcluster labels.  Must exist.
#' @param aggregation \code{"sum"} or \code{"mean"}.  For \code{"sum"},
#'   \code{bundle$counts_raw} is used when available.  Default \code{"mean"}.
#'
#' @return An \code{ObsPointSet} (see \code{\link{ObsPointSet}}).
#' @export
obs_subcluster <- function(bundle, group_col, aggregation = "mean") {
  tryCatch({
    stopifnot(
      is.list(bundle),
      is.matrix(bundle$counts) || inherits(bundle$counts, "Matrix"),
      is.data.frame(bundle$cell_meta),
      is.character(group_col), length(group_col) == 1L,
      group_col %in% names(bundle$cell_meta),
      aggregation %in% c("sum", "mean")
    )

    if (aggregation == "sum") {
      if (is.null(bundle$counts_raw)) {
        warning("obs_subcluster: counts_raw not available; falling back to ",
                "sum of log-normalized counts.")
        count_mat <- bundle$counts
      } else {
        count_mat <- bundle$counts_raw
      }
    } else {
      count_mat <- bundle$counts
    }

    labels    <- as.character(bundle$cell_meta[[group_col]])
    agg_mat   <- as.matrix(.aggregate_cells(count_mat, labels, aggregation))

    lvls      <- colnames(agg_mat)
    n_per     <- setNames(
      as.integer(table(labels)[lvls]),
      lvls
    )
    point_ids <- paste0("subcluster_", lvls)
    colnames(agg_mat) <- point_ids
    names(n_per)      <- point_ids

    point_meta <- data.frame(
      point_id        = point_ids,
      source_subcluster = lvls,
      stringsAsFactors = FALSE
    )

    message("obs_subcluster: ", length(point_ids), " points, median ",
            median(n_per), " cells/point")

    .make_obs(
      matrix      = agg_mat,
      n_cells     = n_per,
      point_meta  = point_meta,
      gene_ids    = rownames(agg_mat),
      design_name = "obs_subcluster",
      design_params = list(group_col = group_col),
      aggregation = aggregation
    )
  }, error = function(e) {
    stop("obs_subcluster failed: ", conditionMessage(e))
  })
}


#' Metacell-style kNN observation points
#'
#' @description
#' Builds metacell-style observation points without requiring external metacell
#' software.  The procedure is:
#' \enumerate{
#'   \item Compute PCA on \code{bundle$counts} (log-normalised).
#'   \item Select \code{n_points} anchor cells via \strong{greedy farthest-point
#'     sampling} in PCA space (see \strong{Caveat} below).
#'   \item For each anchor, collect the anchor itself plus its
#'     \code{target_size - 1} nearest neighbours in PCA space.
#'   \item Aggregate each neighbourhood with \code{\link{.aggregate_cells}}.
#' }
#' Neighbourhoods \strong{may overlap} (\code{overlapping = TRUE} in
#' \code{$design}).
#'
#' @section Caveat:
#' Do \strong{not} use random global subsets as anchors — random bags collapse
#' to the global centroid and produce no spread across cell-state space, which
#' destroys the covariation structure needed for co-expression estimation.  The
#' greedy farthest-point anchor selection is what preserves that structure by
#' ensuring anchors are maximally spread in PCA space.
#'
#' @param bundle      An \code{InputBundle}.
#' @param target_size Integer; approximate number of cells per metacell
#'   (= neighbourhood size, including the anchor).  Default \code{50L}.
#' @param n_points    Integer; number of metacell observation points to create.
#'   Capped to \code{ncol(bundle$counts)} with a warning if larger.  Default
#'   \code{200L}.
#' @param aggregation \code{"sum"} or \code{"mean"}.  For \code{"sum"},
#'   \code{bundle$counts_raw} is used when available.  Default \code{"mean"}.
#'
#' @return An \code{ObsPointSet} (see \code{\link{ObsPointSet}}).
#' @export
obs_metacell_knn <- function(bundle,
                             target_size = 50L,
                             n_points    = 200L,
                             aggregation = "mean") {
  tryCatch({
    stopifnot(
      is.list(bundle),
      is.matrix(bundle$counts) || inherits(bundle$counts, "Matrix"),
      is.data.frame(bundle$cell_meta),
      aggregation %in% c("sum", "mean")
    )
    target_size <- as.integer(target_size)
    n_points    <- as.integer(n_points)
    n_cells     <- ncol(bundle$counts)

    if (n_points > n_cells) {
      warning("obs_metacell_knn: n_points (", n_points, ") > n_cells (",
              n_cells, "); capping to n_cells.")
      n_points <- n_cells
    }

    if (aggregation == "sum") {
      if (is.null(bundle$counts_raw)) {
        warning("obs_metacell_knn: counts_raw not available; falling back to ",
                "sum of log-normalized counts.")
        count_mat <- bundle$counts
      } else {
        count_mat <- bundle$counts_raw
      }
    } else {
      count_mat <- bundle$counts
    }

    # PCA and kNN on log-norm counts
    pca_mat  <- .pca_coords(bundle$counts)
    anchors  <- .greedy_farthest_points(pca_mat, n_points)

    k_nn     <- min(target_size - 1L, n_cells - 1L)

    # Precompute norms for distance calculation
    norms       <- rowSums(pca_mat^2)
    chunk_size  <- min(500L, length(anchors))
    n_a_chunks  <- ceiling(length(anchors) / chunk_size)

    # Build list of neighbour indices per anchor (chunked for memory efficiency)
    nn_lists <- vector("list", length(anchors))
    for (ci in seq_len(n_a_chunks)) {
      idx_start <- (ci - 1L) * chunk_size + 1L
      idx_end   <- min(ci * chunk_size, length(anchors))
      a_idx     <- anchors[idx_start:idx_end]
      batch     <- pca_mat[a_idx, , drop = FALSE]

      sq_dist   <- outer(rowSums(batch^2), norms, "+") - 2 * tcrossprod(batch, pca_mat)
      sq_dist   <- pmax(sq_dist, 0)

      for (li in seq_along(a_idx)) {
        d          <- sq_dist[li, ]
        d[a_idx[li]] <- Inf  # exclude self
        nn         <- order(d)[seq_len(k_nn)]
        nn_lists[[idx_start + li - 1L]] <- c(a_idx[li], nn)
      }
    }

    gene_ids_all <- rownames(bundle$counts)
    n_genes      <- nrow(bundle$counts)
    agg_mat      <- matrix(0.0, nrow = n_genes, ncol = length(anchors))
    rownames(agg_mat) <- gene_ids_all
    point_ids    <- paste0("metacell_", seq_along(anchors))
    colnames(agg_mat) <- point_ids
    n_per        <- integer(length(anchors))

    for (i in seq_along(anchors)) {
      nb    <- nn_lists[[i]]
      sub   <- count_mat[, nb, drop = FALSE]
      n_per[i] <- length(nb)
      if (aggregation == "sum") {
        agg_mat[, i] <- as.numeric(Matrix::rowSums(sub))
      } else {
        agg_mat[, i] <- as.numeric(Matrix::rowMeans(sub))
      }
    }

    n_per_named <- setNames(n_per, point_ids)

    point_meta <- data.frame(
      point_id     = point_ids,
      anchor_index = anchors,
      stringsAsFactors = FALSE
    )

    message("obs_metacell_knn: ", length(anchors), " points, target ",
            target_size, " cells/point")

    .make_obs(
      matrix      = agg_mat,
      n_cells     = n_per_named,
      point_meta  = point_meta,
      gene_ids    = gene_ids_all,
      design_name = "obs_metacell_knn",
      design_params = list(
        target_size = target_size,
        n_points    = length(anchors),
        overlapping = TRUE
      ),
      aggregation = aggregation
    )
  }, error = function(e) {
    stop("obs_metacell_knn failed: ", conditionMessage(e))
  })
}


#' Stratified observation points
#'
#' @description
#' Creates one observation point per unique combination of cell-metadata strata
#' (e.g. cluster × condition).  Combinations with fewer than \code{min_cells}
#' cells are dropped with an informative message.
#'
#' @param bundle      An \code{InputBundle}.
#' @param strata_cols Character vector of column names in \code{bundle$cell_meta}
#'   to cross.  All must exist.
#' @param min_cells   Minimum number of cells a combination must have to be
#'   retained.  Default \code{10L}.
#' @param aggregation \code{"sum"} or \code{"mean"}.  For \code{"sum"},
#'   \code{bundle$counts_raw} is used when available.  Default \code{"mean"}.
#'
#' @return An \code{ObsPointSet} (see \code{\link{ObsPointSet}}).
#' @export
obs_stratified <- function(bundle,
                           strata_cols,
                           min_cells   = 10L,
                           aggregation = "mean") {
  tryCatch({
    stopifnot(
      is.list(bundle),
      is.matrix(bundle$counts) || inherits(bundle$counts, "Matrix"),
      is.data.frame(bundle$cell_meta),
      is.character(strata_cols), length(strata_cols) >= 1L,
      all(strata_cols %in% names(bundle$cell_meta)),
      aggregation %in% c("sum", "mean")
    )
    min_cells <- as.integer(min_cells)

    if (aggregation == "sum") {
      if (is.null(bundle$counts_raw)) {
        warning("obs_stratified: counts_raw not available; falling back to ",
                "sum of log-normalized counts.")
        count_mat <- bundle$counts
      } else {
        count_mat <- bundle$counts_raw
      }
    } else {
      count_mat <- bundle$counts
    }

    cell_meta <- bundle$cell_meta

    # Build composite label
    combo_labels <- apply(
      cell_meta[, strata_cols, drop = FALSE], 1L,
      paste, collapse = "__"
    )

    tbl       <- table(combo_labels)
    n_total   <- length(tbl)
    keep_lvls <- names(tbl)[tbl >= min_cells]
    n_dropped <- n_total - length(keep_lvls)

    if (n_dropped > 0L) {
      message("obs_stratified: dropping ", n_dropped, " strata combination(s) ",
              "with < ", min_cells, " cells.")
    }

    if (length(keep_lvls) == 0L) {
      stop("No strata combinations have >= min_cells (", min_cells, ") cells.")
    }

    # Subset cells to kept combinations
    keep_cells  <- combo_labels %in% keep_lvls
    count_sub   <- count_mat[, keep_cells, drop = FALSE]
    labels_sub  <- combo_labels[keep_cells]

    agg_mat <- as.matrix(.aggregate_cells(count_sub, labels_sub, aggregation))

    lvls      <- colnames(agg_mat)
    n_per     <- setNames(
      as.integer(tbl[lvls]),
      lvls
    )
    point_ids <- paste0("strat_", seq_along(lvls))
    colnames(agg_mat) <- point_ids
    names(n_per)      <- point_ids

    point_meta <- data.frame(
      point_id    = point_ids,
      stratum_key = lvls,
      stringsAsFactors = FALSE
    )
    # Expand stratum key back into individual columns
    split_vals <- strsplit(lvls, "__", fixed = TRUE)
    for (si in seq_along(strata_cols)) {
      point_meta[[strata_cols[si]]] <- vapply(split_vals, `[[`, character(1L), si)
    }

    message("obs_stratified: ", length(point_ids), " points from ",
            n_total, " strata combinations (", n_dropped, " dropped < min_cells)")

    .make_obs(
      matrix      = agg_mat,
      n_cells     = n_per,
      point_meta  = point_meta,
      gene_ids    = rownames(agg_mat),
      design_name = "obs_stratified",
      design_params = list(
        strata_cols = strata_cols,
        min_cells   = min_cells
      ),
      aggregation = aggregation
    )
  }, error = function(e) {
    stop("obs_stratified failed: ", conditionMessage(e))
  })
}


#' Axis-bin observation points
#'
#' @description
#' Creates observation points by binning cells along a continuous axis into
#' equal-frequency (quantile) bins.  The axis can be either:
#' \itemize{
#'   \item A column name already in \code{bundle$cell_meta} holding numeric
#'     values (used directly).
#'   \item A PCA component name such as \code{"PC1"}, \code{"PC2"}, etc. — the
#'     PCA is computed from \code{bundle$counts} and the requested component is
#'     extracted.
#' }
#' Cells with \code{NA} values on the axis are dropped before binning.  Each
#' bin becomes one aggregated observation point.
#'
#' @param bundle      An \code{InputBundle}.
#' @param axis        Character scalar: either a column name in
#'   \code{bundle$cell_meta} or a PCA component (e.g. \code{"PC1"}).
#'   Default \code{"PC1"}.
#' @param n_bins      Integer; number of equal-frequency bins.  Default
#'   \code{20L}.
#' @param aggregation \code{"sum"} or \code{"mean"}.  For \code{"sum"},
#'   \code{bundle$counts_raw} is used when available.  Default \code{"mean"}.
#'
#' @return An \code{ObsPointSet} (see \code{\link{ObsPointSet}}).
#' @export
obs_axis_bin <- function(bundle,
                         axis        = "PC1",
                         n_bins      = 20L,
                         aggregation = "mean") {
  tryCatch({
    stopifnot(
      is.list(bundle),
      is.matrix(bundle$counts) || inherits(bundle$counts, "Matrix"),
      is.data.frame(bundle$cell_meta),
      is.character(axis), length(axis) == 1L,
      aggregation %in% c("sum", "mean")
    )
    n_bins <- as.integer(n_bins)

    if (aggregation == "sum") {
      if (is.null(bundle$counts_raw)) {
        warning("obs_axis_bin: counts_raw not available; falling back to ",
                "sum of log-normalized counts.")
        count_mat <- bundle$counts
      } else {
        count_mat <- bundle$counts_raw
      }
    } else {
      count_mat <- bundle$counts
    }

    cell_meta <- bundle$cell_meta

    # Resolve axis values
    if (axis %in% names(cell_meta)) {
      axis_vals <- as.numeric(cell_meta[[axis]])
    } else {
      # Interpret as PCA component name, e.g. "PC1"
      pc_num <- suppressWarnings(as.integer(sub("^PC", "", axis, ignore.case = TRUE)))
      if (is.na(pc_num) || pc_num < 1L) {
        stop("'axis' ('", axis, "') is neither a column in cell_meta nor a ",
             "valid PCA component name (e.g. 'PC1').")
      }
      pca_mat <- .pca_coords(bundle$counts)
      if (pc_num > ncol(pca_mat)) {
        stop("Requested PC", pc_num, " but PCA only has ", ncol(pca_mat),
             " components.")
      }
      axis_vals <- pca_mat[, pc_num]
    }

    # Drop NAs
    valid     <- !is.na(axis_vals)
    n_dropped <- sum(!valid)
    if (n_dropped > 0L) {
      message("obs_axis_bin: dropping ", n_dropped, " cells with NA on axis '",
              axis, "'.")
    }
    axis_vals_v <- axis_vals[valid]
    count_sub   <- count_mat[, valid, drop = FALSE]

    # Equal-frequency bins via quantile-cut
    probs    <- seq(0, 1, length.out = n_bins + 1L)
    breaks   <- unique(quantile(axis_vals_v, probs = probs))
    if (length(breaks) < 2L) stop("All axis values are identical; cannot bin.")

    bin_labels <- as.character(cut(axis_vals_v, breaks = breaks, include.lowest = TRUE,
                                   labels = FALSE))

    valid_bins  <- !is.na(bin_labels)
    if (!all(valid_bins)) {
      count_sub   <- count_sub[, valid_bins, drop = FALSE]
      bin_labels  <- bin_labels[valid_bins]
    }

    agg_mat <- as.matrix(.aggregate_cells(count_sub, bin_labels, aggregation))

    lvls      <- colnames(agg_mat)
    n_per     <- setNames(
      as.integer(table(bin_labels)[lvls]),
      lvls
    )
    point_ids <- paste0("bin_", lvls)
    colnames(agg_mat) <- point_ids
    names(n_per)      <- point_ids

    # Compute bin midpoint for each level
    lvl_int   <- as.integer(lvls)
    bin_lo    <- breaks[lvl_int]
    bin_hi    <- breaks[lvl_int + 1L]
    bin_mid   <- (bin_lo + bin_hi) / 2

    point_meta <- data.frame(
      point_id    = point_ids,
      bin_index   = as.integer(lvls),
      axis_lo     = bin_lo,
      axis_hi     = bin_hi,
      axis_mid    = bin_mid,
      stringsAsFactors = FALSE
    )

    message("obs_axis_bin: ", length(point_ids), " bins along axis '", axis, "'")

    .make_obs(
      matrix      = agg_mat,
      n_cells     = n_per,
      point_meta  = point_meta,
      gene_ids    = rownames(agg_mat),
      design_name = "obs_axis_bin",
      design_params = list(
        axis   = axis,
        n_bins = n_bins
      ),
      aggregation = aggregation
    )
  }, error = function(e) {
    stop("obs_axis_bin failed: ", conditionMessage(e))
  })
}


# ---------------------------------------------------------------------------
# Normalization and co-expression
# ---------------------------------------------------------------------------

#' Normalize an ObsPointSet matrix
#'
#' @description
#' Apply a normalization method to the \code{$matrix} slot of an
#' \code{\link{ObsPointSet}}.  Returns the \strong{normalized matrix}, not the
#' full ObsPointSet.  Replace the slot with the result:
#' \code{obs$matrix <- normalize_obs(obs)}.
#'
#' @section Depth handling caveat:
#' Spearman correlation is invariant to per-gene monotonic transforms (such as
#' \code{log1p}) but is \strong{not} invariant to per-observation-point depth
#' scaling.  If observation points differ substantially in sequencing depth,
#' depth normalization (\code{"cp10k_log"}) is important for Pearson-based
#' co-expression but may also affect Spearman ranking across points.
#'
#' @param obs    An \code{ObsPointSet}.
#' @param method One of:
#'   \describe{
#'     \item{\code{"none"}}{Return \code{obs$matrix} as-is.}
#'     \item{\code{"cp10k_log"}}{Scale each observation point (column) to 10 000
#'       total counts, then \code{log1p}.  Columns summing to zero are left as
#'       zero.}
#'     \item{\code{"log_only"}}{\code{log1p(obs$matrix)}.}
#'     \item{\code{"zscore_gene"}}{Per-gene z-score across observation points.
#'       Genes with zero standard deviation across points receive a z-score of
#'       zero (not \code{NA} or \code{Inf}).}
#'   }
#'   Default \code{"cp10k_log"}.
#'
#' @return A normalized numeric matrix (genes × observation points).
#' @export
normalize_obs <- function(obs, method = "cp10k_log") {
  stopifnot(
    is.list(obs),
    "matrix" %in% names(obs),
    is.matrix(obs$matrix),
    method %in% c("none", "cp10k_log", "log_only", "zscore_gene")
  )
  mat <- obs$matrix

  switch(method,
    none = mat,

    cp10k_log = {
      col_sums <- colSums(mat)
      scale_f  <- ifelse(col_sums == 0, 1, 1e4 / col_sums)
      scaled   <- sweep(mat, 2L, scale_f, "*")
      log1p(scaled)
    },

    log_only = log1p(mat),

    zscore_gene = {
      gene_means <- rowMeans(mat)
      gene_sds   <- apply(mat, 1L, sd)
      zero_sd    <- gene_sds == 0
      gene_sds[zero_sd] <- 1  # avoid /0; those rows will be 0 after subtract
      z <- sweep(mat - gene_means, 1L, gene_sds, "/")
      z[zero_sd, ] <- 0
      z
    }
  )
}


#' Compute gene-gene co-expression from an ObsPointSet
#'
#' @description
#' Estimates gene-gene co-expression correlations from the aggregated expression
#' matrix in an \code{ObsPointSet}.
#'
#' For Spearman correlation, the rank-transform-then-Pearson approach is used
#' (consistent with \code{estimate_pseudobulk()}): genes are ranked across
#' observation points, and Pearson correlation is computed on those ranks.
#' For Pearson, \code{cor(t(mat))} is computed directly.
#'
#' Only upper-triangle gene pairs with \code{|weight| >= storage_cutoff} are
#' stored in \code{$edge_table} to keep memory manageable.
#'
#' @param obs             An \code{ObsPointSet}.
#' @param cor_type        \code{"spearman"} or \code{"pearson"}.  Default
#'   \code{"spearman"}.
#' @param storage_cutoff  Minimum absolute correlation to retain in the edge
#'   table.  Default \code{0.1}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{$edge_table}}{Data frame: \code{gene_id_A}, \code{gene_id_B},
#'     \code{weight}.  Upper-triangle pairs with \code{|weight| >=
#'     storage_cutoff} only.}
#'   \item{\code{$cor_mat}}{Full genes × genes correlation matrix.}
#'   \item{\code{$gene_ids}}{Character vector of AT-IDs.}
#'   \item{\code{$cor_type}}{The \code{cor_type} argument used.}
#'   \item{\code{$storage_cutoff}}{The \code{storage_cutoff} argument used.}
#' }
#' @export
coexpr_from_obs <- function(obs, cor_type = "spearman", storage_cutoff = 0.1) {
  stopifnot(
    is.list(obs),
    "matrix" %in% names(obs),
    is.matrix(obs$matrix),
    cor_type %in% c("spearman", "pearson"),
    is.numeric(storage_cutoff), length(storage_cutoff) == 1L
  )

  mat      <- obs$matrix
  gene_ids <- obs$gene_ids

  # Spearman via rank-transform-then-Pearson (matches estimate_pseudobulk.R).
  # mat is genes × points (p × N).
  # apply(mat, 1, rank) ranks each gene across points → result is N × p.
  # t(...) restores p × N; cor(t(ranked)) correlates genes (p × p).
  if (cor_type == "spearman") {
    ranked  <- t(apply(mat, 1L, rank))
    cor_mat <- cor(t(ranked))
  } else {
    cor_mat <- cor(t(mat))
  }
  rownames(cor_mat) <- gene_ids
  colnames(cor_mat) <- gene_ids

  # Upper triangle, filtered
  ut     <- which(upper.tri(cor_mat), arr.ind = TRUE)
  w      <- cor_mat[ut]
  keep_e <- abs(w) >= storage_cutoff
  ut     <- ut[keep_e, , drop = FALSE]
  w      <- w[keep_e]

  edge_table <- data.frame(
    gene_id_A = gene_ids[ut[, 1L]],
    gene_id_B = gene_ids[ut[, 2L]],
    weight    = w,
    stringsAsFactors = FALSE
  )

  list(
    edge_table      = edge_table,
    cor_mat         = cor_mat,
    gene_ids        = gene_ids,
    cor_type        = cor_type,
    storage_cutoff  = storage_cutoff
  )
}
