#' @title SingleCellGGM Network Estimation
#'
#' @description
#' Faithful-in-structure reimplementation of Xu, Wang & Ma (2024)
#' Cell Reports Methods 4:100813. Covariance is computed once per stratum;
#' each iteration inverts a 2000-gene submatrix and updates the running
#' minimum |pcor| per pair.
#'
#' @name estimate_singlecellggm
NULL

#' SingleCellGGM co-expression estimation
#'
#' Faithful-in-structure reimplementation of Xu, Wang & Ma (2024)
#' Cell Reports Methods 4:100813 (original MATLAB:
#' github.com/MaShisongLab/SingleCellGGM). Covariance is computed once per
#' stratum; each iteration inverts a 2000-gene submatrix and updates the running
#' minimum |pcor| per pair.
#'
#' @param bundle       InputBundle from load_seurat()
#' @param n_iter       Number of subsampling iterations. NULL (default) =
#'                     auto-compute as round(p*(p-1)/39980) per stratum, so each
#'                     pair is sampled ~100x on average. Only override for testing.
#' @param subsample    Genes per iteration. Default 2000.
#' @param pcor_cutoff  Minimum pcor to retain an edge. Default 0.02.
#' @param coex_cutoff  Minimum number of cells co-detecting both genes. Default 10.
#' @param keep_negative If FALSE (default, matches paper), retain only positive
#'                     partial correlations >= pcor_cutoff. If TRUE, retain edges
#'                     with |pcor| >= pcor_cutoff (signed weight kept).
#' @param ridge        Small value added to covariance submatrix diagonal for
#'                     numerical stability. Default 1e-6.
#' @param seed         Random seed. Default 98.
#' @return Named list of NetworkResult, one per stratum level.
#' @export
estimate_singlecellggm <- function(bundle,
                                   n_iter        = NULL,
                                   subsample     = 2000,
                                   pcor_cutoff   = 0.02,
                                   coex_cutoff   = 10,
                                   keep_negative = FALSE,
                                   ridge         = 1e-6,
                                   seed          = 98) {

  # Warn if reference BLAS detected — performance only, not an error
  blas_lib <- tryCatch(La_library(), error = function(e) "")
  if (nchar(blas_lib) > 0 &&
      !grepl("Accelerate|openblas|mkl|libmkl|blis|atlas|flexiblas",
             blas_lib, ignore.case = TRUE)) {
    message("NOTE: BLAS library (", blas_lib, ") does not appear to be an ",
            "optimised implementation (Accelerate on macOS, OpenBLAS/MKL on Linux). ",
            "The chol2inv() calls will be ~10x slower. Consider relinking R.")
  }

  strat_var <- bundle$stratum_spec$variable
  levels_   <- bundle$stratum_spec$levels
  cell_meta <- bundle$cell_meta
  counts    <- bundle$counts

  results <- list()

  for (level in levels_) {
    idx      <- cell_meta[[strat_var]] == level
    counts_s <- counts[, idx, drop = FALSE]   # genes x cells
    n_cells  <- ncol(counts_s)

    if (n_cells < 2L) {
      warning("Stratum '", level, "': fewer than 2 cells; skipping.")
      next
    }

    # Per-stratum gene filtering:
    #   (a) drop genes detected in < coex_cutoff cells
    n_det    <- Matrix::rowSums(counts_s != 0)
    counts_s <- counts_s[n_det >= coex_cutoff, , drop = FALSE]

    #   (b) drop genes with zero variance (break Cholesky)
    rmu      <- Matrix::rowMeans(counts_s)
    rsq      <- Matrix::rowMeans(counts_s^2)
    keep_var <- (rsq - rmu^2) > .Machine$double.eps
    counts_s <- counts_s[keep_var, , drop = FALSE]

    gene_ids_s <- rownames(counts_s)
    p          <- nrow(counts_s)

    message("SingleCellGGM: stratum '", level, "' — ",
            p, " genes (after per-stratum filter) x ", n_cells, " cells")

    if (p < 2L) {
      warning("Stratum '", level, "': fewer than 2 genes after filtering; skipping.")
      next
    }

    # Resolve effective subsample (cap at p)
    subsample_s <- min(as.integer(subsample), p)

    # Resolve n_iter
    if (is.null(n_iter)) {
      if (subsample_s == p) {
        # Tiny gene set: one full-coverage iteration
        resolved_n_iter <- 1L
        message("  p (", p, ") <= subsample; setting n_iter = 1, subsample = ", p)
      } else {
        resolved_n_iter <- max(1L, round(p * (p - 1L) / 39980))
        message("  auto n_iter = ", resolved_n_iter,
                " (round(", p, "*(", p - 1L, ")/39980); ~100 samplings per pair)")
      }
    } else {
      resolved_n_iter <- as.integer(n_iter)
    }

    # Precompute cov_all and coex ONCE per stratum
    # Densify here; the only time all cells are visited
    counts_dense <- as.matrix(counts_s)
    rm(counts_s)

    # Co-detection matrix: coex[i,j] = # cells detecting both gene i and gene j
    detect <- counts_dense != 0L   # p x n logical
    coex   <- matrix(as.integer(tcrossprod(detect)), nrow = p,
                     dimnames = list(gene_ids_s, gene_ids_s))
    rm(detect)

    # Gene x gene covariance (p x p)
    # Center per gene: subtract per-gene mean. R's column-major recycling of a
    # p-vector over a p x n matrix applies row_means[i] to every element in row i.
    row_means    <- rowMeans(counts_dense)
    centered     <- counts_dense - row_means
    rm(counts_dense)
    cov_all      <- tcrossprod(centered) / (n_cells - 1L)
    rm(centered)
    gc()

    # Initialise accumulators
    pcor_all <- matrix(1.0, nrow = p, ncol = p)  # running min |pcor|, signed
    samp     <- matrix(0L,  nrow = p, ncol = p)  # pair sampling count
    set.seed(seed)

    for (i in seq_len(resolved_n_iter)) {
      if (i %% 1000L == 0L) {
        message("  GGM iteration ", i, " / ", resolved_n_iter)
      }

      j <- sample.int(p, subsample_s)
      S <- cov_all[j, j, drop = FALSE]
      S <- S + diag(ridge, nrow(S))

      # SPD inverse via Cholesky (~2x faster than solve())
      ix <- tryCatch(
        chol2inv(chol(S)),
        error = function(e) {
          # Stronger ridge if Cholesky fails (near-singular block)
          tryCatch(chol2inv(chol(S + diag(1e-4, nrow(S)))),
                   error = function(e2) NULL)
        }
      )
      if (is.null(ix)) next

      dv       <- 1.0 / sqrt(diag(ix))
      pc       <- -(ix * outer(dv, dv))
      diag(pc) <- 1.0

      # Vectorized block min-update — no nested scalar loops
      sub      <- pcor_all[j, j]
      upd      <- abs(pc) < abs(sub)
      sub[upd] <- pc[upd]
      pcor_all[j, j] <- sub
      samp[j, j]     <- samp[j, j] + 1L
    }

    # Zero out pairs never sampled
    pcor_all[samp == 0L] <- 0.0

    # Build edge table from lower triangle (avoid double-counting)
    lt        <- which(lower.tri(pcor_all), arr.ind = TRUE)
    pcor_vals <- pcor_all[lt]
    coex_vals <- coex[lt]
    samp_vals <- samp[lt]

    if (keep_negative) {
      keep_edge <- abs(pcor_vals) >= pcor_cutoff
    } else {
      keep_edge <- pcor_vals >= pcor_cutoff   # positive only — matches paper default
    }
    keep_edge <- keep_edge & (coex_vals >= as.integer(coex_cutoff)) & (samp_vals > 0L)

    if (!any(keep_edge)) {
      warning("Stratum '", level, "': no pairs pass filters (pcor_cutoff=",
              pcor_cutoff, ", coex_cutoff=", coex_cutoff,
              if (keep_negative) ", keep_negative=TRUE" else ", positive-only",
              "); returning empty edge_table.")
      edge_table <- data.frame(
        gene_id_A    = character(0),
        gene_id_B    = character(0),
        weight       = numeric(0),
        coex_cells   = integer(0),
        sampling_num = integer(0),
        stringsAsFactors = FALSE
      )
    } else {
      lt_keep    <- lt[keep_edge, , drop = FALSE]
      # lt columns: [,1] = row (larger matrix index), [,2] = col (smaller)
      # gene_id_A gets the smaller-indexed gene (col), gene_id_B the larger (row)
      edge_table <- data.frame(
        gene_id_A    = gene_ids_s[lt_keep[, 2L]],
        gene_id_B    = gene_ids_s[lt_keep[, 1L]],
        weight       = pcor_vals[keep_edge],
        coex_cells   = as.integer(coex_vals[keep_edge]),
        sampling_num = as.integer(samp_vals[keep_edge]),
        stringsAsFactors = FALSE
      )
    }

    results[[level]] <- list(
      edge_table = edge_table,
      gene_ids   = gene_ids_s,
      stratum_id = level,
      mode       = "singlecellggm",
      params     = list(
        n_iter        = resolved_n_iter,
        subsample     = subsample_s,
        pcor_cutoff   = pcor_cutoff,
        coex_cutoff   = coex_cutoff,
        keep_negative = keep_negative,
        ridge         = ridge,
        seed          = seed,
        aggregation   = "min_abs_pcor_across_iterations",
        n_cells       = n_cells,
        n_genes       = p
      ),
      timestamp = Sys.time()
    )
  }

  results
}
