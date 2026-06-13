#' @title SingleCellGGM Network Estimation
#'
#' @description
#' Estimates a gene co-expression network directly from single-cell counts
#' using the SingleCellGGM method (Xu, Wang & Ma 2024, Cell Reports Methods 4:100813).
#'
#' Algorithm: iterative random subsampling of genes; takes the **minimum** |pcor|
#' across iterations as the conservative final estimate; retains edges with
#' |pcor| >= pcor_cutoff that appeared in >= min_cells iterations.
#'
#' Partial correlation removes indirect edges and partially absorbs the
#' tissue-identity confound by conditioning on other genes. The fundamental
#' paracrine limitation for ligand-receptor pair discovery remains.
#'
#' @name estimate_singlecellggm
NULL

# Fallback partial correlation via MASS::ginv when corpcor is unavailable.
.pcor_ginv <- function(mat) {
  S   <- cov(mat)
  Pi  <- MASS::ginv(S)
  d   <- sqrt(diag(Pi))
  pco <- -Pi / outer(d, d)
  diag(pco) <- 1
  pco
}

#' SingleCellGGM co-expression estimation
#'
#' @param bundle      InputBundle from load_seurat().
#' @param n_iter      Number of subsampling iterations. Default 100.
#' @param subsample   Genes per subsample iteration. Default 2000.
#' @param pcor_cutoff Minimum |pcor| to retain an edge. Default 0.02.
#' @param min_cells   Minimum number of iterations in which a pair must appear.
#'   In the original paper this maps to cells; here it maps to iteration count.
#'   Default 10.
#' @param seed        Random seed. Default 98.
#' @return Named list of NetworkResult, one per stratum level.
#' @export
estimate_singlecellggm <- function(bundle,
                                   n_iter      = 100,
                                   subsample   = 2000,
                                   pcor_cutoff = 0.02,
                                   min_cells   = 10,
                                   seed        = 98) {
  use_corpcor <- requireNamespace("corpcor", quietly = TRUE)
  if (!use_corpcor) {
    warning("Package 'corpcor' not available; using MASS::ginv fallback. ",
            "Install corpcor for shrinkage-based partial correlation.")
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Neither 'corpcor' nor 'MASS' is available. ",
           "Install corpcor (recommended) or MASS.")
    }
  }

  strat_var <- bundle$stratum_spec$variable
  levels_   <- bundle$stratum_spec$levels
  cell_meta <- bundle$cell_meta
  counts    <- bundle$counts

  results <- list()

  for (level in levels_) {
    idx      <- cell_meta[[strat_var]] == level
    counts_s <- counts[, idx, drop = FALSE]
    n_cells  <- ncol(counts_s)
    n_genes  <- nrow(counts_s)
    gene_ids <- rownames(counts_s)
    sub_size <- min(subsample, n_genes)

    message("SingleCellGGM: stratum '", level, "' — ",
            n_genes, " genes x ", n_cells, " cells, ",
            n_iter, " iterations")

    set.seed(seed)

    # Accumulate per-pair observations across iterations.
    # For large sub_size, do.call(rbind, iter_list) can be memory-intensive
    # (sub_size choose 2 rows per iteration). At subsample=2000, n_iter=100,
    # this is ~200M rows; production runs should use data.table or a C++ backend.
    iter_list <- vector("list", n_iter)

    for (i in seq_len(n_iter)) {
      if (i %% 10 == 0) message("  GGM iteration ", i, " / ", n_iter)

      gene_sub_idx <- sample(n_genes, sub_size, replace = FALSE)
      gene_sub     <- gene_ids[gene_sub_idx]
      # Transpose: corpcor/cov expect observations x variables (cells x genes)
      mat_sub      <- t(as.matrix(counts_s[gene_sub_idx, , drop = FALSE]))

      pcor_mat <- if (use_corpcor) {
        tryCatch(
          as.matrix(corpcor::pcor.shrink(mat_sub, verbose = FALSE)),
          error = function(e) {
            warning("corpcor::pcor.shrink failed in iteration ", i,
                    ": ", conditionMessage(e), ". Using MASS::ginv fallback.")
            .pcor_ginv(mat_sub)
          }
        )
      } else {
        .pcor_ginv(mat_sub)
      }

      ut <- which(upper.tri(pcor_mat), arr.ind = TRUE)
      iter_list[[i]] <- data.frame(
        ga   = gene_sub[ut[, 1L]],
        gb   = gene_sub[ut[, 2L]],
        pcor = pcor_mat[ut],
        stringsAsFactors = FALSE
      )
    }

    all_pairs <- do.call(rbind, iter_list)

    if (is.null(all_pairs) || nrow(all_pairs) == 0L) {
      warning("Stratum '", level, "': no pairs accumulated; skipping.")
      next
    }

    # Canonical ordering: ensure gene_id_A < gene_id_B lexicographically
    swap <- all_pairs$ga > all_pairs$gb
    if (any(swap)) {
      tmp               <- all_pairs$ga[swap]
      all_pairs$ga[swap] <- all_pairs$gb[swap]
      all_pairs$gb[swap] <- tmp
    }
    all_pairs$abs_pcor <- abs(all_pairs$pcor)
    all_pairs$key      <- paste(all_pairs$ga, all_pairs$gb, sep = "\x1f")

    # Aggregation: min |pcor| across iterations, signed by mean pcor direction
    agg_n       <- tapply(all_pairs$abs_pcor, all_pairs$key, length)
    agg_sum     <- tapply(all_pairs$pcor,     all_pairs$key, sum)
    agg_min_abs <- tapply(all_pairs$abs_pcor, all_pairs$key, min)

    pair_names   <- names(agg_n)
    n_app        <- as.integer(agg_n)
    sum_pcor     <- as.numeric(agg_sum)
    min_abs_pcor <- as.numeric(agg_min_abs)
    final_weight <- sign(sum_pcor) * min_abs_pcor

    # Filters: |pcor| >= pcor_cutoff AND appeared in >= min_cells iterations
    keep <- (min_abs_pcor >= pcor_cutoff) & (n_app >= min_cells)

    if (!any(keep)) {
      warning("Stratum '", level, "': no pairs pass filters (pcor_cutoff=",
              pcor_cutoff, ", min_cells=", min_cells, "); returning empty edge_table.")
      edge_table <- data.frame(
        gene_id_A = character(0),
        gene_id_B = character(0),
        weight    = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      kept_names <- pair_names[keep]
      split_keys <- strsplit(kept_names, "\x1f", fixed = TRUE)
      edge_table <- data.frame(
        gene_id_A = vapply(split_keys, `[[`, character(1L), 1L),
        gene_id_B = vapply(split_keys, `[[`, character(1L), 2L),
        weight    = final_weight[keep],
        stringsAsFactors = FALSE
      )
    }

    results[[level]] <- list(
      edge_table = edge_table,
      gene_ids   = gene_ids,
      stratum_id = level,
      mode       = "singlecellggm",
      params     = list(
        n_iter      = n_iter,
        subsample   = subsample,
        pcor_cutoff = pcor_cutoff,
        min_cells   = min_cells,
        seed        = seed,
        aggregation = "min_abs_pcor_across_iterations",
        n_cells     = n_cells,
        n_genes     = n_genes
      ),
      timestamp = Sys.time()
    )
  }

  results
}
