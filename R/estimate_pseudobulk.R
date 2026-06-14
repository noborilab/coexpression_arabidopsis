#' @title Pseudobulk Network Estimation
#'
#' @description
#' Aggregates cells to pseudobulk replicates per stratum, then estimates a
#' marginal Spearman correlation network (rank-transform-then-Pearson).
#' Produces one NetworkResult per stratum level.
#'
#' Robustness statistics (R_score) are computed downstream in `R/robustness.R`.
#'
#' @name estimate_pseudobulk
NULL

#' Pseudobulk co-expression estimation
#'
#' @param bundle      InputBundle from load_seurat().
#' @param min_samples Minimum number of pseudobulk samples required per stratum
#'   to run correlation. Default 5.
#' @param min_expr    Minimum mean expression (in pseudobulk) for a gene to be
#'   retained per stratum. Default 0 (keep all after adapter filter).
#' @return Named list of NetworkResult, one per stratum level.
#' @export
estimate_pseudobulk <- function(bundle, min_samples = 5, min_expr = 0) {
  strat_var <- bundle$stratum_spec$variable
  levels_   <- bundle$stratum_spec$levels
  cell_meta <- bundle$cell_meta
  counts    <- bundle$counts

  grp_col <- if ("group_var" %in% names(cell_meta)) "group_var" else NULL
  if (is.null(grp_col)) {
    warning("No 'group_var' column in cell_meta. ",
            "Treating each cell as its own pseudobulk unit. ",
            "This is unusual and may produce unreliable correlations.")
  }

  # Storage filter: practical limit on edge table size (not a significance filter)
  storage_cutoff <- 0.1

  # First pass: aggregate cells to pseudobulk per stratum
  n_cells_per_stratum      <- setNames(integer(length(levels_)), levels_)
  n_pseudobulk_per_stratum <- setNames(integer(length(levels_)), levels_)
  pseudobulk_list          <- vector("list", length(levels_))
  names(pseudobulk_list)   <- levels_

  for (level in levels_) {
    idx <- cell_meta[[strat_var]] == level
    n_cells_per_stratum[level] <- sum(idx)

    if (!any(idx)) next

    counts_s <- counts[, idx, drop = FALSE]
    meta_s   <- cell_meta[idx, , drop = FALSE]

    if (!is.null(grp_col)) {
      groups        <- as.character(meta_s[[grp_col]])
      unique_groups <- unique(groups)

      pb <- vapply(unique_groups, function(g) {
        ci <- groups == g
        if (sum(ci) == 1L) {
          as.numeric(counts_s[, ci, drop = TRUE])
        } else {
          as.numeric(Matrix::rowMeans(counts_s[, ci, drop = FALSE]))
        }
      }, numeric(nrow(counts_s)))

      rownames(pb) <- rownames(counts_s)
      colnames(pb) <- unique_groups
    } else {
      pb <- as.matrix(counts_s)
    }

    n_pseudobulk_per_stratum[level] <- ncol(pb)
    pseudobulk_list[[level]]        <- pb
  }

  # Second pass: filter, correlate, assemble NetworkResults
  results <- list()

  for (level in levels_) {
    pb <- pseudobulk_list[[level]]

    if (is.null(pb)) {
      warning("Stratum '", level, "': no cells found; skipping.")
      next
    }

    if (ncol(pb) < min_samples) {
      warning("Stratum '", level, "': ", ncol(pb),
              " pseudobulk samples < min_samples = ", min_samples, "; skipping.")
      next
    }

    keep <- rowMeans(pb) >= min_expr
    pb   <- pb[keep, , drop = FALSE]

    if (nrow(pb) == 0L) {
      warning("Stratum '", level, "': no genes pass min_expr filter; skipping.")
      next
    }

    gene_ids_used <- rownames(pb)

    # Spearman via rank-transform-then-Pearson (see docs/BACKGROUND.md).
    # pb is genes × samples (p × N).
    # apply(pb, 1, rank) ranks each gene ACROSS samples → result is N × p.
    # t(...)  restores p × N so each row is one gene's rank profile across samples.
    # cor(t(ranked)) correlates columns of the N × p matrix = gene-gene Spearman (p × p).
    ranked  <- t(apply(pb, 1, rank))   # p × N: gene i ranked across N samples
    cor_mat <- cor(t(ranked))

    # Upper triangle only, storage filter |weight| >= 0.1
    ut     <- which(upper.tri(cor_mat), arr.ind = TRUE)
    w      <- cor_mat[ut]
    keep_e <- abs(w) >= storage_cutoff
    ut     <- ut[keep_e, , drop = FALSE]
    w      <- w[keep_e]

    edge_table <- data.frame(
      gene_id_A = gene_ids_used[ut[, 1L]],
      gene_id_B = gene_ids_used[ut[, 2L]],
      weight    = w,
      stringsAsFactors = FALSE
    )

    results[[level]] <- list(
      edge_table = edge_table,
      gene_ids   = gene_ids_used,
      stratum_id = level,
      mode       = "pseudobulk",
      params     = list(
        min_samples              = min_samples,
        min_expr                 = min_expr,
        group_var                = grp_col,
        storage_cutoff           = storage_cutoff,
        n_cells_per_stratum      = n_cells_per_stratum,
        n_pseudobulk_per_stratum = n_pseudobulk_per_stratum
      ),
      timestamp = Sys.time()
    )
  }

  results
}
