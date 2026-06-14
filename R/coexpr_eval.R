#' @title Prior-Free Observation-Point Design Evaluation
#'
#' @description
#' Diagnostic metrics for evaluating and comparing observation-point designs
#' without reference to gold-standard gene sets. All metrics are computable
#' from a single ObsPointSet (or from an InputBundle + design function).
#'
#' **Metric overview:**
#' - `eval_effective_rank` — participation ratio of singular values; measures
#'   how many independent covariation axes the design resolves.
#' - `eval_visible_genes` — fraction of genes with non-degenerate variance.
#' - `eval_heldout_predictivity` — cross-validated guilt-by-association R²;
#'   the single most informative prior-free metric.
#' - `eval_null_gap` — fraction of suprathreshold edges vs. permuted null.
#' - `eval_depth_leakage` — correlation of gene network degree with mean
#'   expression; detects abundance confounding.
#' - `eval_splithalf` — split-half reproducibility of the correlation matrix.
#' - `eval_downsample_depth` — robustness to binomial depth thinning.
#' - `eval_downsample_cells` — robustness to cell-count reduction.
#' - `evaluate_obs_design` — full harness combining all of the above.
#'
#' @name coexpr_eval
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Subset an InputBundle to a subset of columns (cells).
# cell_indices: integer vector of column positions in bundle$counts.
.subset_bundle_cells <- function(bundle, cell_indices) {
  new_bundle <- bundle
  new_bundle$counts    <- bundle$counts[, cell_indices, drop = FALSE]
  if (!is.null(bundle$counts_raw)) {
    new_bundle$counts_raw <- bundle$counts_raw[, cell_indices, drop = FALSE]
  }
  new_bundle$cell_meta <- bundle$cell_meta[cell_indices, , drop = FALSE]
  new_bundle
}

# Binomial thinning of a raw-count matrix.
# fraction in (0, 1]: each count c -> rbinom(1, c, fraction).
# Returns an integer matrix of the same dimensions, or NULL if counts_raw is NULL.
.thin_counts <- function(counts_raw, fraction) {
  if (is.null(counts_raw)) return(NULL)
  if (fraction >= 1.0)     return(counts_raw)
  mat_int <- as.matrix(counts_raw)
  thinned <- matrix(
    rbinom(length(mat_int), as.integer(mat_int), prob = fraction),
    nrow = nrow(mat_int), ncol = ncol(mat_int),
    dimnames = dimnames(mat_int)
  )
  storage.mode(thinned) <- "integer"
  thinned
}

# Build an ObsPointSet from a bundle using a design function, then normalize.
# design_fn is called as design_fn(bundle, design_args[[1]], ...) via do.call.
# normalize_obs(obs, norm_method) replaces obs$matrix with the normalized version.
# Wraps design_fn call in tryCatch to surface design failures clearly.
.build_obs_from_bundle <- function(bundle, design_fn, design_args, norm_method) {
  obs <- tryCatch(
    do.call(design_fn, c(list(bundle), design_args)),
    error = function(e) stop("Design function failed: ", conditionMessage(e))
  )
  obs$matrix <- normalize_obs(obs, norm_method)
  obs
}

# Compute genes×genes correlation matrix from a gene×points matrix.
# Spearman: rank-transform each gene across points, then Pearson on t(ranked).
# Pearson: direct cor(t(mat)).
# No storage filter — returns the full square matrix.
.cor_mat_from_obs <- function(obs, cor_type) {
  mat <- obs$matrix
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    m <- matrix(NA_real_, nrow = nrow(mat), ncol = nrow(mat))
    rownames(m) <- colnames(m) <- rownames(mat)
    return(m)
  }
  if (cor_type == "spearman") {
    ranked <- t(apply(mat, 1L, rank))
    cor(t(ranked))
  } else {
    cor(t(mat))
  }
}

# Jaccard similarity of top-k% edges (by |cor|) between two correlation matrices.
# k = max(1, round(n_upper_triangle * k_frac)).
# Returns scalar in [0, 1].
.jaccard_topk <- function(m1, m2, k_frac = 0.01) {
  ut      <- upper.tri(m1, diag = FALSE)
  n_upper <- sum(ut)
  k       <- max(1L, round(n_upper * k_frac))

  v1 <- m1[ut]
  v2 <- m2[ut]

  top1 <- order(abs(v1), decreasing = TRUE)[seq_len(min(k, n_upper))]
  top2 <- order(abs(v2), decreasing = TRUE)[seq_len(min(k, n_upper))]

  n_intersect <- length(intersect(top1, top2))
  n_union     <- length(union(top1, top2))

  if (n_union == 0L) return(NA_real_)
  n_intersect / n_union
}

# Pearson correlation of the vectorized upper triangles of two square matrices.
# Returns NA if either vector has zero variance.
.matrix_correlation <- function(m1, m2) {
  ut <- upper.tri(m1, diag = FALSE)
  v1 <- m1[ut]
  v2 <- m2[ut]
  if (stats::var(v1) == 0 || stats::var(v2) == 0) return(NA_real_)
  stats::cor(v1, v2)
}

# Per-gene network degree: number of edges with |cor| >= threshold.
# Returns a numeric vector of length nrow(cor_mat).
.gene_network_degree <- function(cor_mat, threshold = 0.3) {
  diag(cor_mat) <- 0
  rowSums(abs(cor_mat) >= threshold)
}

# ---------------------------------------------------------------------------
# Exported eval functions
# ---------------------------------------------------------------------------

#' Split-half reproducibility of a co-expression design
#'
#' @description
#' Randomly split cells in `bundle` into two equal halves `n_reps` times. For
#' each split, independently build an ObsPointSet and correlation matrix in each
#' half (using `.build_obs_from_bundle`). Agreement is measured by Pearson
#' correlation of the upper-triangle elements of the two correlation matrices
#' (`mat_cor`) and Jaccard of the top-`top_k_frac` edges (`jaccard`). Only
#' genes common to both halves are compared.
#'
#' If a split fails (e.g. too few cells for the design), that rep is skipped
#' with a warning and the effective `n_reps` is reduced. If no rep succeeds,
#' all metric columns are `NA`.
#'
#' @param bundle InputBundle (see package description).
#' @param design_fn Function that maps a bundle to an ObsPointSet.
#' @param design_args Named list of extra arguments forwarded to `design_fn`
#'   after the bundle. Default `list()`.
#' @param cor_type Correlation type: `"spearman"` (default) or `"pearson"`.
#' @param norm_method Normalization method passed to `normalize_obs()`.
#'   Default `"cp10k_log"`.
#' @param n_reps Number of random split-half replicates. Default `5L`.
#' @param top_k_frac Fraction of edges used for Jaccard comparison. Default `0.01`.
#' @return A one-row `data.frame` with columns: `mat_cor_mean`, `mat_cor_sd`,
#'   `jaccard_mean`, `jaccard_sd`, `n_reps` (effective number of successful reps).
#' @export
eval_splithalf <- function(bundle,
                           design_fn,
                           design_args    = list(),
                           cor_type       = "spearman",
                           norm_method    = "cp10k_log",
                           n_reps         = 5L,
                           top_k_frac     = 0.01) {

  n_cells    <- ncol(bundle$counts)
  mat_cors   <- numeric(0L)
  jaccards   <- numeric(0L)

  for (rep in seq_len(n_reps)) {
    idx      <- sample.int(n_cells)
    half     <- floor(n_cells / 2L)
    idx1     <- idx[seq_len(half)]
    idx2     <- idx[seq(half + 1L, n_cells)]

    obs1 <- tryCatch(
      .build_obs_from_bundle(.subset_bundle_cells(bundle, idx1),
                             design_fn, design_args, norm_method),
      error = function(e) {
        warning("Split-half rep ", rep, " half-1 failed: ", conditionMessage(e))
        NULL
      }
    )
    obs2 <- tryCatch(
      .build_obs_from_bundle(.subset_bundle_cells(bundle, idx2),
                             design_fn, design_args, norm_method),
      error = function(e) {
        warning("Split-half rep ", rep, " half-2 failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(obs1) || is.null(obs2)) next

    common_genes <- intersect(obs1$gene_ids, obs2$gene_ids)
    if (length(common_genes) < 2L) {
      warning("Split-half rep ", rep, ": fewer than 2 common genes; skipping.")
      next
    }

    cm1 <- .cor_mat_from_obs(
      structure(list(matrix = obs1$matrix[common_genes, , drop = FALSE]),
                class = "list"),
      cor_type
    )
    cm2 <- .cor_mat_from_obs(
      structure(list(matrix = obs2$matrix[common_genes, , drop = FALSE]),
                class = "list"),
      cor_type
    )

    mc <- tryCatch(.matrix_correlation(cm1, cm2), error = function(e) NA_real_)
    jc <- tryCatch(.jaccard_topk(cm1, cm2, k_frac = top_k_frac),
                   error = function(e) NA_real_)

    mat_cors <- c(mat_cors, mc)
    jaccards <- c(jaccards, jc)
  }

  n_ok <- length(mat_cors)

  if (n_ok == 0L) {
    warning("eval_splithalf: no successful reps; returning NA.")
    return(data.frame(
      mat_cor_mean  = NA_real_,
      mat_cor_sd    = NA_real_,
      jaccard_mean  = NA_real_,
      jaccard_sd    = NA_real_,
      n_reps        = 0L
    ))
  }

  data.frame(
    mat_cor_mean  = mean(mat_cors, na.rm = TRUE),
    mat_cor_sd    = if (n_ok > 1L) stats::sd(mat_cors, na.rm = TRUE) else NA_real_,
    jaccard_mean  = mean(jaccards, na.rm = TRUE),
    jaccard_sd    = if (n_ok > 1L) stats::sd(jaccards, na.rm = TRUE) else NA_real_,
    n_reps        = n_ok
  )
}

# ---------------------------------------------------------------------------

#' Depth-downsampling robustness
#'
#' @description
#' Binomially thin `bundle$counts_raw` to each fraction in `fractions`,
#' rebuild the ObsPointSet and correlation matrix, and compare each result to
#' the full-depth result (fraction = 1.0 if included in `fractions`, otherwise
#' the result from the unthinned bundle).
#'
#' **Note on thinning and normalization:** This function thins `counts_raw`
#' only. When the design uses `aggregation = "sum"`, the thinned `counts_raw`
#' propagates correctly through the sum-aggregation path. When the design uses
#' `aggregation = "mean"`, the log-normalized `counts` matrix is unchanged, so
#' thinning has no effect; a warning is emitted in this case if detected via
#' `design_args$aggregation`.
#'
#' @param bundle InputBundle. Must have non-NULL `counts_raw`; if `NULL`,
#'   the function warns and returns `NULL`.
#' @param design_fn Function that maps a bundle to an ObsPointSet.
#' @param design_args Named list of extra arguments forwarded to `design_fn`.
#'   Default `list()`.
#' @param fractions Numeric vector of thinning fractions in (0, 1].
#'   Default `c(1.0, 0.5, 0.25, 0.1)`.
#' @param cor_type Correlation type: `"spearman"` (default) or `"pearson"`.
#' @param norm_method Normalization method passed to `normalize_obs()`.
#'   Default `"cp10k_log"`.
#' @return A `data.frame` with columns `fraction`, `mat_cor`, `jaccard`,
#'   or `NULL` if `bundle$counts_raw` is `NULL`.
#' @export
eval_downsample_depth <- function(bundle,
                                  design_fn,
                                  design_args = list(),
                                  fractions   = c(1.0, 0.5, 0.25, 0.1),
                                  cor_type    = "spearman",
                                  norm_method = "cp10k_log") {

  if (is.null(bundle$counts_raw)) {
    warning("eval_downsample_depth: bundle$counts_raw is NULL. ",
            "Depth downsampling requires raw integer counts. Returning NULL.")
    return(NULL)
  }

  agg <- design_args[["aggregation"]]
  if (!is.null(agg) && identical(agg, "mean")) {
    warning("eval_downsample_depth: design_args$aggregation = 'mean'. ",
            "The log-normalized counts matrix is not thinned, so depth ",
            "thinning will have no effect for mean-aggregation designs.")
  }

  fractions <- sort(unique(fractions), decreasing = TRUE)

  # Build the reference at fraction = 1.0 (full depth).
  ref_obs <- tryCatch(
    .build_obs_from_bundle(bundle, design_fn, design_args, norm_method),
    error = function(e) {
      warning("eval_downsample_depth: full-depth build failed: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(ref_obs)) {
    return(data.frame(fraction = fractions, mat_cor = NA_real_,
                      jaccard  = NA_real_))
  }
  ref_cm <- .cor_mat_from_obs(ref_obs, cor_type)

  out_rows <- vector("list", length(fractions))

  for (fi in seq_along(fractions)) {
    frac <- fractions[fi]

    if (frac >= 1.0) {
      out_rows[[fi]] <- data.frame(fraction = frac,
                                   mat_cor  = 1.0,
                                   jaccard  = 1.0)
      next
    }

    thinned_raw <- .thin_counts(bundle$counts_raw, frac)
    mod_bundle  <- bundle
    mod_bundle$counts_raw <- thinned_raw

    obs_t <- tryCatch(
      .build_obs_from_bundle(mod_bundle, design_fn, design_args, norm_method),
      error = function(e) {
        warning("eval_downsample_depth: fraction ", frac, " failed: ",
                conditionMessage(e))
        NULL
      }
    )

    if (is.null(obs_t)) {
      out_rows[[fi]] <- data.frame(fraction = frac,
                                   mat_cor  = NA_real_,
                                   jaccard  = NA_real_)
      next
    }

    common_genes <- intersect(ref_obs$gene_ids, obs_t$gene_ids)

    if (length(common_genes) < 2L) {
      out_rows[[fi]] <- data.frame(fraction = frac,
                                   mat_cor  = NA_real_,
                                   jaccard  = NA_real_)
      next
    }

    cm_t <- .cor_mat_from_obs(
      structure(list(matrix = obs_t$matrix[common_genes, , drop = FALSE]),
                class = "list"),
      cor_type
    )
    cm_r <- ref_cm[common_genes, common_genes, drop = FALSE]

    mc <- tryCatch(.matrix_correlation(cm_r, cm_t), error = function(e) NA_real_)
    jc <- tryCatch(.jaccard_topk(cm_r, cm_t), error = function(e) NA_real_)

    out_rows[[fi]] <- data.frame(fraction = frac, mat_cor = mc, jaccard = jc)
  }

  do.call(rbind, out_rows)
}

# ---------------------------------------------------------------------------

#' Cell-downsampling robustness
#'
#' @description
#' Sample cells to each fraction in `fractions`, rebuild the ObsPointSet and
#' correlation matrix, and compare each result to the full-cell result
#' (fraction = 1.0 if included in `fractions`, otherwise the result from the
#' full bundle). Agreement is measured by Pearson correlation of upper-triangle
#' elements and Jaccard of the top-1% edges.
#'
#' @param bundle InputBundle (see package description).
#' @param design_fn Function that maps a bundle to an ObsPointSet.
#' @param design_args Named list of extra arguments forwarded to `design_fn`.
#'   Default `list()`.
#' @param fractions Numeric vector of cell-fraction levels in (0, 1].
#'   Default `c(1.0, 0.75, 0.5, 0.25)`.
#' @param cor_type Correlation type: `"spearman"` (default) or `"pearson"`.
#' @param norm_method Normalization method passed to `normalize_obs()`.
#'   Default `"cp10k_log"`.
#' @return A `data.frame` with columns `fraction`, `mat_cor`, `jaccard`.
#' @export
eval_downsample_cells <- function(bundle,
                                  design_fn,
                                  design_args = list(),
                                  fractions   = c(1.0, 0.75, 0.5, 0.25),
                                  cor_type    = "spearman",
                                  norm_method = "cp10k_log") {

  fractions <- sort(unique(fractions), decreasing = TRUE)
  n_cells   <- ncol(bundle$counts)

  # Build full-depth reference.
  ref_obs <- tryCatch(
    .build_obs_from_bundle(bundle, design_fn, design_args, norm_method),
    error = function(e) {
      warning("eval_downsample_cells: full-cell build failed: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(ref_obs)) {
    return(data.frame(fraction = fractions, mat_cor = NA_real_,
                      jaccard  = NA_real_))
  }
  ref_cm <- .cor_mat_from_obs(ref_obs, cor_type)

  out_rows <- vector("list", length(fractions))

  for (fi in seq_along(fractions)) {
    frac <- fractions[fi]

    if (frac >= 1.0) {
      out_rows[[fi]] <- data.frame(fraction = frac,
                                   mat_cor  = 1.0,
                                   jaccard  = 1.0)
      next
    }

    n_keep    <- max(1L, round(n_cells * frac))
    keep_idx  <- sample.int(n_cells, size = n_keep, replace = FALSE)
    sub_bundle <- .subset_bundle_cells(bundle, keep_idx)

    obs_s <- tryCatch(
      .build_obs_from_bundle(sub_bundle, design_fn, design_args, norm_method),
      error = function(e) {
        warning("eval_downsample_cells: fraction ", frac, " failed: ",
                conditionMessage(e))
        NULL
      }
    )

    if (is.null(obs_s)) {
      out_rows[[fi]] <- data.frame(fraction = frac,
                                   mat_cor  = NA_real_,
                                   jaccard  = NA_real_)
      next
    }

    common_genes <- intersect(ref_obs$gene_ids, obs_s$gene_ids)

    if (length(common_genes) < 2L) {
      out_rows[[fi]] <- data.frame(fraction = frac,
                                   mat_cor  = NA_real_,
                                   jaccard  = NA_real_)
      next
    }

    cm_s <- .cor_mat_from_obs(
      structure(list(matrix = obs_s$matrix[common_genes, , drop = FALSE]),
                class = "list"),
      cor_type
    )
    cm_r <- ref_cm[common_genes, common_genes, drop = FALSE]

    mc <- tryCatch(.matrix_correlation(cm_r, cm_s), error = function(e) NA_real_)
    jc <- tryCatch(.jaccard_topk(cm_r, cm_s), error = function(e) NA_real_)

    out_rows[[fi]] <- data.frame(fraction = frac, mat_cor = mc, jaccard = jc)
  }

  do.call(rbind, out_rows)
}

# ---------------------------------------------------------------------------

#' Effective rank of an observation-point matrix
#'
#' @description
#' Computes the participation ratio of singular values of the gene-centered
#' observation-point matrix: `eff_rank = (sum(s_i))^2 / sum(s_i^2)`. This is
#' a prior-free measure of the number of independent covariation axes that the
#' observation points resolve. Higher effective rank indicates that the design
#' captures a broader range of coexpression structure.
#'
#' If the matrix has rank 0 (all-zero after centering), `eff_rank` is returned
#' as 0.
#'
#' @param obs ObsPointSet (see package description). Uses `obs$matrix`
#'   (genes × observation-points).
#' @return A one-row `data.frame` with columns: `eff_rank` (numeric),
#'   `n_points` (integer), `n_genes` (integer).
#' @export
eval_effective_rank <- function(obs) {
  mat <- obs$matrix

  n_genes  <- nrow(mat)
  n_points <- ncol(mat)

  if (n_genes == 0L || n_points == 0L) {
    return(data.frame(eff_rank = 0, n_points = n_points, n_genes = n_genes))
  }

  # Center each gene across observation points.
  mat_centered <- mat - rowMeans(mat)

  sv <- tryCatch(
    svd(mat_centered, nu = 0L, nv = 0L)$d,
    error = function(e) numeric(0L)
  )

  sv_pos <- sv[sv > 0]

  if (length(sv_pos) == 0L) {
    return(data.frame(eff_rank = 0, n_points = n_points, n_genes = n_genes))
  }

  eff_rank <- (sum(sv_pos))^2 / sum(sv_pos^2)

  data.frame(
    eff_rank = eff_rank,
    n_points = as.integer(n_points),
    n_genes  = as.integer(n_genes)
  )
}

# ---------------------------------------------------------------------------

#' Count genes with non-degenerate variance across observation points
#'
#' @description
#' A gene is "visible" if its variance across observation points exceeds
#' `min_var` AND it is not all-zero. Invisible genes contribute no information
#' to co-expression estimates and should be inspected as a design quality
#' indicator.
#'
#' @param obs ObsPointSet (see package description).
#' @param min_var Minimum variance threshold. Default `1e-6`.
#' @return A one-row `data.frame` with columns: `n_visible` (integer),
#'   `n_total` (integer), `frac_visible` (numeric).
#' @export
eval_visible_genes <- function(obs, min_var = 1e-6) {
  mat     <- obs$matrix
  n_total <- nrow(mat)

  if (n_total == 0L || ncol(mat) == 0L) {
    return(data.frame(n_visible   = 0L,
                      n_total     = n_total,
                      frac_visible = NA_real_))
  }

  gene_vars    <- apply(mat, 1L, stats::var)
  gene_nonzero <- rowSums(mat != 0) > 0L

  n_visible <- sum(gene_vars > min_var & gene_nonzero, na.rm = TRUE)

  data.frame(
    n_visible    = as.integer(n_visible),
    n_total      = as.integer(n_total),
    frac_visible = n_visible / n_total
  )
}

# ---------------------------------------------------------------------------

#' Cross-validated guilt-by-association predictivity
#'
#' @description
#' Partition observation points into `n_folds` folds. For each fold, compute
#' gene-gene correlations on training points, identify the top-`k_partners`
#' correlated partners for each gene (excluding self), predict held-out
#' expression as the absolute-value-weighted mean of partner expression, and
#' compute per-gene R². Returns mean and median R² across genes.
#'
#' This is the single most informative prior-free metric: it rewards both real
#' AND broad covariation structure simultaneously.
#'
#' **Implementation notes:**
#' - `set.seed` is not called internally; callers may set the seed before
#'   calling.
#' - If `n_points < 2 * n_folds`, `n_folds` is reduced to
#'   `max(2L, floor(n_points / 3L))`.
#' - Only genes with `ss_tot > 0` in the held-out fold contribute to R².
#' - Prediction: `pred_g = sum(|w_gj| * expr_j) / sum(|w_gj|)` for the top-k
#'   partners j.
#'
#' @param obs ObsPointSet (see package description).
#' @param k_partners Number of top-correlated partner genes to use for
#'   prediction. Default `10L`.
#' @param n_folds Number of cross-validation folds. Default `5L`.
#' @param cor_type Correlation type for computing partner correlations:
#'   `"spearman"` (default) or `"pearson"`.
#' @return A one-row `data.frame` with columns: `predictivity_mean_r2`,
#'   `predictivity_median_r2`.
#' @export
eval_heldout_predictivity <- function(obs,
                                      k_partners = 10L,
                                      n_folds    = 5L,
                                      cor_type   = "spearman") {

  mat      <- obs$matrix
  n_genes  <- nrow(mat)
  n_points <- ncol(mat)

  na_result <- data.frame(predictivity_mean_r2   = NA_real_,
                          predictivity_median_r2  = NA_real_)

  if (n_genes < 2L || n_points < 4L) {
    warning("eval_heldout_predictivity: too few genes or points; returning NA.")
    return(na_result)
  }

  # Large-n guard: building an n_genes×n_genes weight matrix is O(n²) in memory
  # and causes GC blowup above ~500 genes. Subsample genes for the estimate.
  max_genes_full <- 500L
  gene_idx <- seq_len(n_genes)
  if (n_genes > max_genes_full) {
    set.seed(42L)
    gene_idx <- sort(sample(n_genes, max_genes_full))
    mat      <- mat[gene_idx, , drop = FALSE]
    n_genes  <- max_genes_full
    message("eval_heldout_predictivity: subsampling to ", max_genes_full,
            " genes for O(n²) guard (full n_genes = ", nrow(obs$matrix), ").")
  }

  # Reduce n_folds if too few points.
  if (n_points < 2L * n_folds) {
    n_folds <- max(2L, floor(n_points / 3L))
  }

  fold_ids <- sample(rep(seq_len(n_folds), length.out = n_points))
  r2_mat   <- matrix(NA_real_, nrow = n_genes, ncol = n_folds)

  k_use <- min(k_partners, n_genes - 1L)

  for (fold in seq_len(n_folds)) {
    test_idx  <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)

    if (length(train_idx) < 2L || length(test_idx) < 1L) next

    mat_train <- mat[, train_idx, drop = FALSE]
    mat_test  <- mat[, test_idx,  drop = FALSE]

    cor_train <- .cor_mat_from_obs(
      structure(list(matrix = mat_train), class = "list"),
      cor_type
    )
    diag(cor_train) <- 0

    # Compute predictions gene-by-gene to avoid materialising an n²-element
    # weight matrix. Each gene's top-k partners are found from one row of
    # cor_train, and the prediction is a weighted mean of partner values.
    n_test <- length(test_idx)
    preds  <- matrix(NA_real_, nrow = n_genes, ncol = n_test)

    for (g in seq_len(n_genes)) {
      row_abs <- abs(cor_train[g, ])
      top_k_idx <- order(row_abs, decreasing = TRUE)[seq_len(k_use)]
      w <- row_abs[top_k_idx]
      w_sum <- sum(w)
      if (w_sum == 0) next
      preds[g, ] <- as.numeric((w / w_sum) %*% mat_test[top_k_idx, , drop = FALSE])
    }

    rm(cor_train); gc(verbose = FALSE)

    # Per-gene R² on held-out fold; clip to [-1, 1] to guard against
    # numerical blow-up from near-zero ss_tot.
    y_mean <- rowMeans(mat_test)
    ss_res <- rowSums((mat_test - preds)^2, na.rm = TRUE)
    ss_tot <- rowSums((mat_test - y_mean)^2)

    r2_fold <- ifelse(ss_tot == 0, NA_real_, 1 - ss_res / ss_tot)
    r2_fold <- pmax(-1, pmin(1, r2_fold))   # clip; removes blow-up outliers
    r2_mat[, fold] <- r2_fold
  }

  # Average R² across folds per gene, then summarise across genes.
  gene_r2 <- rowMeans(r2_mat, na.rm = TRUE)
  gene_r2 <- gene_r2[is.finite(gene_r2)]

  if (length(gene_r2) == 0L) return(na_result)

  data.frame(
    predictivity_mean_r2   = mean(gene_r2),
    predictivity_median_r2 = stats::median(gene_r2)
  )
}

# ---------------------------------------------------------------------------

#' Real vs. permuted-null edge-weight distribution
#'
#' @description
#' Compare the fraction of gene-pair edges with `|cor| >= threshold` in real
#' data vs. an independently-shuffled null. For each of `n_perm` permutations,
#' the rows of `obs$matrix` are independently permuted (destroying covariation
#' structure), correlations are recomputed, and the suprathreshold fraction is
#' recorded. The null-gap ratio (`real_frac / perm_frac_mean`) is the primary
#' output: values well above 1 indicate genuine coexpression signal.
#'
#' @param obs ObsPointSet (see package description).
#' @param cor_type Correlation type: `"spearman"` (default) or `"pearson"`.
#' @param n_perm Number of permutations for the null distribution. Default `20L`.
#' @param threshold Absolute-correlation threshold for edge counting.
#'   Default `0.3`.
#' @return A one-row `data.frame` with columns: `null_gap_ratio`
#'   (real_frac / perm_frac_mean; `Inf` if perm_frac_mean = 0),
#'   `real_frac`, `perm_frac_mean`.
#' @export
eval_null_gap <- function(obs,
                          cor_type  = "spearman",
                          n_perm    = 20L,
                          threshold = 0.3) {

  mat <- obs$matrix

  na_result <- data.frame(null_gap_ratio = NA_real_,
                          real_frac      = NA_real_,
                          perm_frac_mean = NA_real_)

  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    warning("eval_null_gap: matrix has fewer than 2 genes or 2 points; ",
            "returning NA.")
    return(na_result)
  }

  .count_frac <- function(m, cor_tp) {
    cm <- .cor_mat_from_obs(structure(list(matrix = m), class = "list"), cor_tp)
    ut <- upper.tri(cm, diag = FALSE)
    mean(abs(cm[ut]) >= threshold, na.rm = TRUE)
  }

  real_frac <- .count_frac(mat, cor_type)

  perm_fracs <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    # apply(mat, 1, sample) returns a points×genes matrix; transpose to genes×points.
    mat_perm      <- t(apply(mat, 1L, sample))
    perm_fracs[p] <- .count_frac(mat_perm, cor_type)
  }

  perm_frac_mean  <- mean(perm_fracs, na.rm = TRUE)
  null_gap_ratio  <- if (perm_frac_mean == 0) Inf else real_frac / perm_frac_mean

  data.frame(
    null_gap_ratio = null_gap_ratio,
    real_frac      = real_frac,
    perm_frac_mean = perm_frac_mean
  )
}

# ---------------------------------------------------------------------------

#' Depth / abundance confounding diagnostic
#'
#' @description
#' Computes the Spearman correlation between each gene's network degree (number
#' of edges with `|cor| >= threshold`) and its mean expression across
#' observation points. A high positive correlation indicates that the design or
#' normalization is leaking a sequencing-depth / abundance axis: hubs are
#' simply highly expressed genes rather than genuinely co-regulated ones.
#'
#' @param obs ObsPointSet (see package description).
#' @param threshold Absolute-correlation threshold used to define network edges.
#'   Default `0.3`.
#' @return A one-row `data.frame` with columns: `depth_leakage_rho` (Spearman
#'   rho between degree and mean expression), `n_genes` (integer).
#' @export
eval_depth_leakage <- function(obs, threshold = 0.3) {
  mat     <- obs$matrix
  n_genes <- nrow(mat)

  na_result <- data.frame(depth_leakage_rho = NA_real_,
                          n_genes           = as.integer(n_genes))

  if (n_genes < 3L || ncol(mat) < 2L) {
    warning("eval_depth_leakage: too few genes or points; returning NA.")
    return(na_result)
  }

  cm     <- .cor_mat_from_obs(obs, cor_type = "spearman")
  degree <- .gene_network_degree(cm, threshold = threshold)
  mean_expr <- rowMeans(mat)

  rho <- tryCatch(
    stats::cor(degree, mean_expr, method = "spearman"),
    error = function(e) NA_real_
  )

  data.frame(
    depth_leakage_rho = rho,
    n_genes           = as.integer(n_genes)
  )
}

# ---------------------------------------------------------------------------

#' Full prior-free evaluation harness for one observation-point design
#'
#' @description
#' **IMPORTANT CAVEAT:** Stability metrics (split-half, downsampling) alone
#' favour trivial designs (one huge observation point, or designs that only
#' resolve the dominant axis). They **MUST** be read jointly with richness
#' metrics (effective rank, visible genes, held-out predictivity). Selection
#' should be made on the stability–richness Pareto front, not any single
#' metric.
#'
#' Runs the full suite of prior-free evaluation metrics for a given
#' observation-point design:
#' - Richness: `eval_effective_rank`, `eval_visible_genes`,
#'   `eval_heldout_predictivity`
#' - Signal: `eval_null_gap`, `eval_depth_leakage`
#' - Stability: `eval_splithalf` (controlled by `n_splithalf` / `splithalf_reps`)
#' - Optional: `eval_downsample_depth`, `eval_downsample_cells`
#'
#' Each metric is wrapped in `tryCatch`; failures record `NA` with a warning
#' so that one failing metric does not abort the entire evaluation.
#'
#' @param bundle InputBundle (see package description).
#' @param design_fn Function that maps a bundle to an ObsPointSet.
#' @param design_args Named list of extra arguments forwarded to `design_fn`.
#'   Default `list()`.
#' @param cor_type Correlation type: `"spearman"` (default) or `"pearson"`.
#' @param norm_method Normalization method passed to `normalize_obs()`.
#'   Default `"cp10k_log"`.
#' @param n_splithalf Number of split-half replicates. Pass `0L` to skip
#'   split-half evaluation. Default `3L`.
#' @param run_downsample_depth Logical; if `TRUE`, run `eval_downsample_depth`
#'   and append `downsample_depth_mat_cor_at_0.5`. Default `FALSE`.
#' @param run_downsample_cells Logical; if `TRUE`, run `eval_downsample_cells`
#'   and append `downsample_cells_mat_cor_at_0.5`. Default `FALSE`.
#' @param splithalf_reps Alias for `n_splithalf` (kept for API compatibility).
#'   If both are provided, `n_splithalf` takes precedence. Default `3L`.
#' @param heldout_folds Number of cross-validation folds for
#'   `eval_heldout_predictivity`. Default `5L`.
#' @param null_perm Number of permutations for `eval_null_gap`. Default `10L`.
#' @return A one-row `data.frame` with columns: `design_name`, `n_points`,
#'   `n_genes`, `eff_rank`, `n_visible`, `frac_visible`,
#'   `predictivity_mean_r2`, `predictivity_median_r2`, `null_gap_ratio`,
#'   `real_frac`, `perm_frac_mean`, `depth_leakage_rho`,
#'   `splithalf_mat_cor_mean`, `splithalf_mat_cor_sd`,
#'   `splithalf_jaccard_mean`, `splithalf_jaccard_sd`.
#'   Plus `downsample_depth_mat_cor_at_0.5` if `run_downsample_depth = TRUE`.
#'   Plus `downsample_cells_mat_cor_at_0.5` if `run_downsample_cells = TRUE`.
#' @export
evaluate_obs_design <- function(bundle,
                                design_fn,
                                design_args           = list(),
                                cor_type              = "spearman",
                                norm_method           = "cp10k_log",
                                n_splithalf           = 3L,
                                run_downsample_depth  = FALSE,
                                run_downsample_cells  = FALSE,
                                splithalf_reps        = 3L,
                                heldout_folds         = 5L,
                                null_perm             = 10L) {

  # splithalf_reps is an alias; n_splithalf takes precedence.
  effective_splithalf <- n_splithalf

  # ---- Build the full ObsPointSet once for cheap metrics ------------------

  obs <- tryCatch(
    .build_obs_from_bundle(bundle, design_fn, design_args, norm_method),
    error = function(e) {
      warning("evaluate_obs_design: design build failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(obs)) {
    warning("evaluate_obs_design: could not build ObsPointSet; ",
            "returning all-NA row.")
    empty <- data.frame(
      design_name            = NA_character_,
      n_points               = NA_integer_,
      n_genes                = NA_integer_,
      eff_rank               = NA_real_,
      n_visible              = NA_integer_,
      frac_visible           = NA_real_,
      predictivity_mean_r2   = NA_real_,
      predictivity_median_r2 = NA_real_,
      null_gap_ratio         = NA_real_,
      real_frac              = NA_real_,
      perm_frac_mean         = NA_real_,
      depth_leakage_rho      = NA_real_,
      splithalf_mat_cor_mean = NA_real_,
      splithalf_mat_cor_sd   = NA_real_,
      splithalf_jaccard_mean = NA_real_,
      splithalf_jaccard_sd   = NA_real_,
      stringsAsFactors = FALSE
    )
    if (run_downsample_depth)
      empty$downsample_depth_mat_cor_at_0.5 <- NA_real_
    if (run_downsample_cells)
      empty$downsample_cells_mat_cor_at_0.5 <- NA_real_
    return(empty)
  }

  design_name <- obs$design$name %||% NA_character_
  n_points    <- as.integer(ncol(obs$matrix))
  n_genes     <- as.integer(nrow(obs$matrix))

  # ---- Effective rank -----------------------------------------------------

  eff_rank_res <- tryCatch(
    eval_effective_rank(obs),
    error = function(e) {
      warning("evaluate_obs_design [eff_rank]: ", conditionMessage(e))
      data.frame(eff_rank = NA_real_, n_points = n_points, n_genes = n_genes)
    }
  )

  # ---- Visible genes -------------------------------------------------------

  vis_res <- tryCatch(
    eval_visible_genes(obs),
    error = function(e) {
      warning("evaluate_obs_design [visible_genes]: ", conditionMessage(e))
      data.frame(n_visible = NA_integer_, n_total = n_genes,
                 frac_visible = NA_real_)
    }
  )

  # ---- Held-out predictivity -----------------------------------------------

  pred_res <- tryCatch(
    eval_heldout_predictivity(obs, n_folds = heldout_folds,
                              cor_type = cor_type),
    error = function(e) {
      warning("evaluate_obs_design [heldout_predictivity]: ",
              conditionMessage(e))
      data.frame(predictivity_mean_r2   = NA_real_,
                 predictivity_median_r2 = NA_real_)
    }
  )

  # ---- Null gap ------------------------------------------------------------

  null_res <- tryCatch(
    eval_null_gap(obs, cor_type = cor_type, n_perm = null_perm),
    error = function(e) {
      warning("evaluate_obs_design [null_gap]: ", conditionMessage(e))
      data.frame(null_gap_ratio = NA_real_, real_frac = NA_real_,
                 perm_frac_mean = NA_real_)
    }
  )

  # ---- Depth leakage -------------------------------------------------------

  leak_res <- tryCatch(
    eval_depth_leakage(obs),
    error = function(e) {
      warning("evaluate_obs_design [depth_leakage]: ", conditionMessage(e))
      data.frame(depth_leakage_rho = NA_real_, n_genes = n_genes)
    }
  )

  # ---- Split-half ----------------------------------------------------------

  if (effective_splithalf > 0L) {
    sh_res <- tryCatch(
      eval_splithalf(bundle,
                     design_fn   = design_fn,
                     design_args = design_args,
                     cor_type    = cor_type,
                     norm_method = norm_method,
                     n_reps      = effective_splithalf),
      error = function(e) {
        warning("evaluate_obs_design [splithalf]: ", conditionMessage(e))
        data.frame(mat_cor_mean = NA_real_, mat_cor_sd = NA_real_,
                   jaccard_mean = NA_real_, jaccard_sd = NA_real_,
                   n_reps       = 0L)
      }
    )
  } else {
    sh_res <- data.frame(mat_cor_mean = NA_real_, mat_cor_sd = NA_real_,
                         jaccard_mean = NA_real_, jaccard_sd = NA_real_,
                         n_reps       = 0L)
  }

  # ---- Assemble output row -------------------------------------------------

  out <- data.frame(
    design_name            = design_name,
    n_points               = n_points,
    n_genes                = n_genes,
    eff_rank               = eff_rank_res$eff_rank,
    n_visible              = vis_res$n_visible,
    frac_visible           = vis_res$frac_visible,
    predictivity_mean_r2   = pred_res$predictivity_mean_r2,
    predictivity_median_r2 = pred_res$predictivity_median_r2,
    null_gap_ratio         = null_res$null_gap_ratio,
    real_frac              = null_res$real_frac,
    perm_frac_mean         = null_res$perm_frac_mean,
    depth_leakage_rho      = leak_res$depth_leakage_rho,
    splithalf_mat_cor_mean = sh_res$mat_cor_mean,
    splithalf_mat_cor_sd   = sh_res$mat_cor_sd,
    splithalf_jaccard_mean = sh_res$jaccard_mean,
    splithalf_jaccard_sd   = sh_res$jaccard_sd,
    stringsAsFactors = FALSE
  )

  # ---- Optional: downsample depth ------------------------------------------

  if (run_downsample_depth) {
    dd_res <- tryCatch(
      eval_downsample_depth(bundle,
                            design_fn   = design_fn,
                            design_args = design_args,
                            fractions   = c(1.0, 0.5, 0.25, 0.1),
                            cor_type    = cor_type,
                            norm_method = norm_method),
      error = function(e) {
        warning("evaluate_obs_design [downsample_depth]: ", conditionMessage(e))
        NULL
      }
    )
    val_at_half <- NA_real_
    if (!is.null(dd_res)) {
      row_half <- dd_res[abs(dd_res$fraction - 0.5) < 1e-9, , drop = FALSE]
      if (nrow(row_half) > 0L) val_at_half <- row_half$mat_cor[1L]
    }
    out$downsample_depth_mat_cor_at_0.5 <- val_at_half
  }

  # ---- Optional: downsample cells ------------------------------------------

  if (run_downsample_cells) {
    dc_res <- tryCatch(
      eval_downsample_cells(bundle,
                            design_fn   = design_fn,
                            design_args = design_args,
                            fractions   = c(1.0, 0.75, 0.5, 0.25),
                            cor_type    = cor_type,
                            norm_method = norm_method),
      error = function(e) {
        warning("evaluate_obs_design [downsample_cells]: ", conditionMessage(e))
        NULL
      }
    )
    val_at_half <- NA_real_
    if (!is.null(dc_res)) {
      row_half <- dc_res[abs(dc_res$fraction - 0.5) < 1e-9, , drop = FALSE]
      if (nrow(row_half) > 0L) val_at_half <- row_half$mat_cor[1L]
    }
    out$downsample_cells_mat_cor_at_0.5 <- val_at_half
  }

  out
}
