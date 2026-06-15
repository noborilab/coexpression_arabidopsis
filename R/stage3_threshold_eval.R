#' @title Stage 3 Threshold Evaluation Adapters (obs-point basis)
#'
#' @description
#' Prior-free evaluation of thresholded co-expression networks using the SAME
#' obs-point axis as Stage 1/2 (genes × 298 pseudobulk profiles), not the
#' 4-condition fingerprint matrix used in the prior (invalid) Phase 2.
#'
#' All five functions share the same input contract:
#' - `obs`:      ObsPointSet produced by obs_subcluster + normalize_obs.
#'               `obs$matrix` must be genes × n_obs_points (genes in rownames).
#' - `edges_dt`: data.table with at minimum columns `gene_id_A`, `gene_id_B`,
#'               and `mean_abs_r` (the retaining threshold's basis).
#' - `min_abs_r` / `top_k`: exactly one non-NULL (the threshold lever).
#'
#' **No data.table GForce operations anywhere** — this machine (aarch64-darwin)
#' segfaults on frank(by=), dt[group,.()], unique(by=), and dt[bool_expr] on
#' large tables. All operations use base-R: rank(), apply(), ave(), which(),
#' rowSums(), rowMeans(), svd(), cor(), and plain matrix indexing.
#'
#' @name stage3_threshold_eval
NULL

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

#' Compute a genes×genes Spearman correlation matrix from a column subset.
#'
#' Uses BLAS matrix multiplication on rank-normalized rows, avoiding the O(n²)
#' memory overhead of computing all gene-pair correlations one by one.
#' Input: `mat` (genes × n_pts), `col_idx` (integer column positions to use).
#' Returns: symmetric genes×genes matrix of Spearman correlations.
#' @keywords internal
.s3_spearman_mat <- function(mat, col_idx) {
  sub <- mat[, col_idx, drop = FALSE]
  # Rank each gene (row) independently across the selected obs-points
  rk  <- t(apply(sub, 1L, rank))      # genes × n_cols
  # Center rows and normalize to unit length
  rc  <- rk - rowMeans(rk)
  ss  <- sqrt(rowSums(rc^2))
  ss[ss < 1e-12] <- 1e-12
  rn  <- rc / ss                       # unit-row matrix
  # Spearman ≡ Pearson of ranked = rn %*% t(rn) via BLAS
  rn %*% t(rn)
}

#' Extract pair correlations from a square matrix using integer index pairs.
#' @keywords internal
.s3_extract_pairs <- function(cor_mat, iA, iB) cor_mat[cbind(iA, iB)]

#' Per-gene top-k union mask (base-R ave, no data.table frank/GForce).
#' Returns logical: TRUE if pair is in top-k by |r| for gene_A OR gene_B.
#' @keywords internal
.s3_topk_mask <- function(abs_r_vals, gA, gB, k) {
  rk_A <- ave(-abs_r_vals, gA, FUN = rank)
  rk_B <- ave(-abs_r_vals, gB, FUN = rank)
  ((rk_A <= k) | (rk_B <= k)) & !is.na(abs_r_vals)
}

#' Subset obs$matrix to visible genes and build integer indices for candidate pairs.
#' Returns list: mat, gA, gB, iA, iB (only pairs where both genes are in obs).
#' @keywords internal
.s3_net_subset <- function(obs, edges_dt) {
  mat_full <- obs$matrix
  net_g    <- sort(unique(c(edges_dt$gene_id_A, edges_dt$gene_id_B)))
  avail    <- rownames(mat_full) %in% net_g
  mat      <- mat_full[avail, , drop = FALSE]
  g_names  <- rownames(mat)
  gidx     <- setNames(seq_len(nrow(mat)), g_names)

  ok <- edges_dt$gene_id_A %in% g_names & edges_dt$gene_id_B %in% g_names
  gA <- edges_dt$gene_id_A[ok]
  gB <- edges_dt$gene_id_B[ok]
  iA <- gidx[gA]
  iB <- gidx[gB]

  list(mat = mat, gA = gA, gB = gB, iA = iA, iB = iB,
       n_missing = sum(!ok))
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. stage3_eval_splithalf
# ─────────────────────────────────────────────────────────────────────────────

#' Split-half stability on the obs-point axis
#'
#' For each rep: randomly splits the n_obs_points obs-point columns into two
#' equal halves, computes the genes×genes Spearman correlation matrix from each
#' half independently, applies the threshold (global |r| or per-gene top-k),
#' and reports Jaccard of the retained edge sets (primary) plus Pearson of the
#' full upper-triangle correlation structures (secondary).
#'
#' @param obs ObsPointSet; obs$matrix must be genes × n_obs_points with gene
#'   IDs as rownames.
#' @param edges_dt data.table: gene_id_A, gene_id_B, mean_abs_r.
#' @param min_abs_r Numeric; global |r| threshold (Lever A). Exactly one of
#'   min_abs_r or top_k must be non-NULL.
#' @param top_k Integer; per-gene top-k threshold (Lever B).
#' @param n_reps Integer; number of random split-half replicates. Default 5.
#' @param seed Integer RNG seed. Default 98.
#' @return data.frame: splithalf_jaccard, splithalf_pearson, splithalf_jaccard_sd,
#'   n_reps (effective reps completed).
#' @export
stage3_eval_splithalf <- function(obs, edges_dt,
                                  min_abs_r = NULL, top_k    = NULL,
                                  n_reps    = 5L,   seed     = 98L) {
  if (is.null(min_abs_r) == is.null(top_k))
    stop("stage3_eval_splithalf: exactly one of min_abs_r or top_k must be set.")

  ns    <- .s3_net_subset(obs, edges_dt)
  mat   <- ns$mat
  gA    <- ns$gA;  gB <- ns$gB
  iA    <- ns$iA;  iB <- ns$iB
  n_pts <- ncol(mat)
  half  <- floor(n_pts / 2L)
  n_genes <- nrow(mat)

  if (nrow(mat) < 2L || n_pts < 4L || length(iA) == 0L)
    return(data.frame(splithalf_jaccard    = NA_real_,
                      splithalf_pearson    = NA_real_,
                      splithalf_jaccard_sd = NA_real_,
                      n_reps               = 0L))

  set.seed(seed)
  jacs  <- numeric(n_reps)
  pears <- numeric(n_reps)
  n_ok  <- 0L

  for (rep in seq_len(n_reps)) {
    perm  <- sample.int(n_pts)
    cols1 <- perm[seq_len(half)]
    cols2 <- perm[seq(half + 1L, n_pts)]

    cm1 <- .s3_spearman_mat(mat, cols1)
    cm2 <- .s3_spearman_mat(mat, cols2)

    r1 <- .s3_extract_pairs(cm1, iA, iB)
    r2 <- .s3_extract_pairs(cm2, iA, iB)

    if (!is.null(min_abs_r)) {
      in1 <- !is.na(r1) & abs(r1) >= min_abs_r
      in2 <- !is.na(r2) & abs(r2) >= min_abs_r
    } else {
      in1 <- .s3_topk_mask(abs(r1), gA, gB, top_k)
      in2 <- .s3_topk_mask(abs(r2), gA, gB, top_k)
    }

    n_int <- sum(in1 & in2)
    n_uni <- sum(in1 | in2)
    jac   <- if (n_uni == 0L) NA_real_ else n_int / n_uni

    # Pearson of upper-tri: proxy via candidate pairs (memory-safe for large n)
    # For n_genes <= 5000 use full upper-tri; otherwise use candidate-pair proxy.
    if (n_genes <= 5000L) {
      ut   <- upper.tri(cm1, diag = FALSE)
      pear <- tryCatch(stats::cor(cm1[ut], cm2[ut]),
                       error = function(e) NA_real_)
    } else {
      pear <- tryCatch(stats::cor(r1, r2), error = function(e) NA_real_)
    }

    rm(cm1, cm2); gc(verbose = FALSE)

    if (!is.na(jac)) {
      n_ok        <- n_ok + 1L
      jacs[n_ok]  <- jac
      pears[n_ok] <- if (is.finite(pear)) pear else NA_real_
    }
  }

  if (n_ok == 0L)
    return(data.frame(splithalf_jaccard    = NA_real_,
                      splithalf_pearson    = NA_real_,
                      splithalf_jaccard_sd = NA_real_,
                      n_reps               = 0L))

  data.frame(
    splithalf_jaccard    = mean(jacs[seq_len(n_ok)]),
    splithalf_pearson    = mean(pears[seq_len(n_ok)], na.rm = TRUE),
    splithalf_jaccard_sd = if (n_ok > 1L) stats::sd(jacs[seq_len(n_ok)]) else NA_real_,
    n_reps               = as.integer(n_ok)
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. stage3_eval_effective_rank
# ─────────────────────────────────────────────────────────────────────────────

#' Effective rank of the obs-point matrix masked to visible genes
#'
#' Subsets obs$matrix to genes with ≥1 retained edge at this threshold and
#' computes the participation ratio of singular values of the masked genes×298
#' matrix: eff_rank = (Σ sᵢ)² / Σ sᵢ².  Can range up to min(n_visible, 298).
#'
#' @param obs ObsPointSet.
#' @param edges_dt data.table: gene_id_A, gene_id_B, mean_abs_r.
#' @return data.frame: eff_rank, n_visible, n_points.
#' @export
stage3_eval_effective_rank <- function(obs, edges_dt) {
  net_genes <- unique(c(edges_dt$gene_id_A, edges_dt$gene_id_B))
  mat_full  <- obs$matrix
  avail     <- rownames(mat_full) %in% net_genes
  mat_vis   <- mat_full[avail, , drop = FALSE]

  n_vis <- nrow(mat_vis)
  n_pts <- ncol(mat_vis)

  if (n_vis < 2L || n_pts < 2L)
    return(data.frame(eff_rank  = NA_real_,
                      n_visible = as.integer(n_vis),
                      n_points  = as.integer(n_pts)))

  mat_c <- mat_vis - rowMeans(mat_vis)    # center genes across obs-points
  sv    <- tryCatch(svd(mat_c, nu = 0L, nv = 0L)$d,
                    error = function(e) numeric(0L))
  sv    <- sv[sv > 0]

  er <- if (length(sv) == 0L) NA_real_ else (sum(sv))^2 / sum(sv^2)

  data.frame(eff_rank  = er,
             n_visible = as.integer(n_vis),
             n_points  = as.integer(n_pts))
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. stage3_eval_heldout
# ─────────────────────────────────────────────────────────────────────────────

#' Cross-validated guilt-by-association R² on the obs-point axis
#'
#' 5-fold CV partitioning the 298 obs-point columns.  For each fold:
#' (1) compute Spearman on training obs-points, (2) apply the threshold to get
#' retained edges, (3) for each gene predict its test-fold expression from its
#' retained network neighbours via weighted mean (weights = |train Spearman|),
#' (4) compute per-gene R² clipped to [-1, 1].  Returns the mean across all
#' folds and all genes with at least one retained neighbour.
#'
#' @param obs ObsPointSet.
#' @param edges_dt data.table: gene_id_A, gene_id_B, mean_abs_r.
#' @param min_abs_r Numeric or NULL (Lever A).
#' @param top_k Integer or NULL (Lever B).
#' @param n_folds Integer; CV folds. Default 5.
#' @param seed Integer RNG seed. Default 98.
#' @return data.frame: heldout_r2.
#' @export
stage3_eval_heldout <- function(obs, edges_dt,
                                min_abs_r = NULL, top_k    = NULL,
                                n_folds   = 5L,   seed     = 98L) {
  if (is.null(min_abs_r) == is.null(top_k))
    stop("stage3_eval_heldout: exactly one of min_abs_r or top_k must be set.")

  ns      <- .s3_net_subset(obs, edges_dt)
  mat     <- ns$mat
  gA      <- ns$gA;  gB <- ns$gB
  iA      <- ns$iA;  iB <- ns$iB
  n_genes <- nrow(mat);  n_pts <- ncol(mat)

  if (n_genes < 2L || n_pts < 4L || length(iA) == 0L)
    return(data.frame(heldout_r2 = NA_real_))

  set.seed(seed)
  fold_ids <- sample(rep(seq_len(n_folds), length.out = n_pts))
  all_r2   <- numeric(0L)

  for (fold in seq_len(n_folds)) {
    train_idx <- which(fold_ids != fold)
    test_idx  <- which(fold_ids == fold)
    if (length(train_idx) < 2L || length(test_idx) < 1L) next

    # Spearman on training columns
    cm_tr   <- .s3_spearman_mat(mat, train_idx)
    r_tr    <- .s3_extract_pairs(cm_tr, iA, iB)
    rm(cm_tr); gc(verbose = FALSE)

    # Apply threshold
    if (!is.null(min_abs_r)) {
      retained <- !is.na(r_tr) & abs(r_tr) >= min_abs_r
    } else {
      retained <- .s3_topk_mask(abs(r_tr), gA, gB, top_k)
    }
    if (!any(retained, na.rm = TRUE)) next

    iA_r <- iA[retained];  iB_r <- iB[retained]
    w_r  <- abs(r_tr[retained])

    mat_test <- mat[, test_idx, drop = FALSE]
    r2_fold  <- rep(NA_real_, n_genes)

    # Build adjacency: for gene g as target, collect (source, weight) via split()
    # Both directions: predict B from A, and A from B.
    tgt_all <- c(iB_r, iA_r)
    src_all <- c(iA_r, iB_r)
    w_all   <- c(w_r,  w_r)

    # Group by target gene (base-R split — no data.table)
    grps <- split(
      data.frame(src = src_all, w = w_all, stringsAsFactors = FALSE),
      tgt_all
    )

    for (g_str in names(grps)) {
      g      <- as.integer(g_str)
      rec    <- grps[[g_str]]
      w_sum  <- sum(rec$w)
      if (w_sum < 1e-10) next

      w_norm <- rec$w / w_sum
      pred   <- as.numeric(w_norm %*% mat_test[rec$src, , drop = FALSE])
      actual <- mat_test[g, ]
      ss_tot <- sum((actual - mean(actual))^2)
      if (ss_tot < 1e-10) next

      r2_fold[g] <- pmax(-1, pmin(1, 1 - sum((actual - pred)^2) / ss_tot))
    }

    rm(mat_test, grps); gc(verbose = FALSE)
    all_r2 <- c(all_r2, r2_fold[is.finite(r2_fold)])
  }

  r2 <- if (length(all_r2) > 0L) mean(all_r2, na.rm = TRUE) else NA_real_
  if (!is.finite(r2)) r2 <- NA_real_
  data.frame(heldout_r2 = r2)
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. stage3_eval_null_gap
# ─────────────────────────────────────────────────────────────────────────────

#' Real vs. permuted-null edge density on the obs-point axis
#'
#' Computes the Spearman correlation matrix from the full obs-point set for the
#' visible gene subset, counts the fraction of candidate pairs exceeding the
#' threshold (real_frac), then repeats on independently row-permuted obs matrices
#' (destroying cross-gene correlations while preserving per-gene distributions).
#' null_gap = real_frac / mean(perm_fracs); values >> 1 indicate genuine signal.
#'
#' @param obs ObsPointSet.
#' @param edges_dt data.table: gene_id_A, gene_id_B, mean_abs_r.
#' @param min_abs_r Numeric or NULL.
#' @param top_k Integer or NULL.
#' @param n_perm Integer; permutations. Default 10.
#' @param seed Integer. Default 98.
#' @return data.frame: null_gap, real_frac, perm_frac_mean.
#' @export
stage3_eval_null_gap <- function(obs, edges_dt,
                                 min_abs_r = NULL, top_k    = NULL,
                                 n_perm    = 10L,  seed     = 98L) {
  if (is.null(min_abs_r) == is.null(top_k))
    stop("stage3_eval_null_gap: exactly one of min_abs_r or top_k must be set.")

  ns     <- .s3_net_subset(obs, edges_dt)
  mat    <- ns$mat
  gA     <- ns$gA;  gB <- ns$gB
  iA     <- ns$iA;  iB <- ns$iB
  n_pts  <- ncol(mat)

  if (nrow(mat) < 2L || n_pts < 2L || length(iA) == 0L)
    return(data.frame(null_gap       = NA_real_,
                      real_frac      = NA_real_,
                      perm_frac_mean = NA_real_))

  .count_frac <- function(m) {
    cm <- .s3_spearman_mat(m, seq_len(ncol(m)))
    r  <- .s3_extract_pairs(cm, iA, iB)
    rm(cm)
    if (!is.null(min_abs_r)) {
      mean(abs(r) >= min_abs_r, na.rm = TRUE)
    } else {
      mean(.s3_topk_mask(abs(r), gA, gB, top_k), na.rm = TRUE)
    }
  }

  real_frac <- .count_frac(mat)

  set.seed(seed)
  perm_fracs <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    # Independently shuffle each gene's expression across obs-points (base-R apply)
    mat_perm   <- t(apply(mat, 1L, sample))
    perm_fracs[p] <- .count_frac(mat_perm)
    rm(mat_perm); gc(verbose = FALSE)
  }

  pf_mean  <- mean(perm_fracs)
  null_gap <- if (pf_mean < 1e-10) Inf else real_frac / pf_mean

  data.frame(null_gap       = null_gap,
             real_frac      = real_frac,
             perm_frac_mean = pf_mean)
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. stage3_eval_visible_genes
# ─────────────────────────────────────────────────────────────────────────────

#' Non-isolated gene count for a thresholded network
#'
#' @param edges_dt data.table with gene_id_A and gene_id_B.
#' @param n_total Integer; total gene universe size. Default 11010.
#' @return data.frame: n_visible, n_total, frac_visible.
#' @export
stage3_eval_visible_genes <- function(edges_dt, n_total = 11010L) {
  if (nrow(edges_dt) == 0L)
    return(data.frame(n_visible    = 0L,
                      n_total      = as.integer(n_total),
                      frac_visible = 0))

  n_vis <- length(unique(c(edges_dt$gene_id_A, edges_dt$gene_id_B)))
  data.frame(n_visible    = as.integer(n_vis),
             n_total      = as.integer(n_total),
             frac_visible = n_vis / n_total)
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. stage3_eval_louvain (descriptive only — NOT a selection metric)
# ─────────────────────────────────────────────────────────────────────────────

#' Fast Louvain module statistics (descriptive support, not selection)
#'
#' @param edges_dt data.table: gene_id_A, gene_id_B, mean_abs_r (weight).
#' @param seed Integer. Default 98.
#' @return data.frame: n_modules, grey_rate, median_module_size.
#' @export
stage3_eval_louvain <- function(edges_dt, seed = 98L) {
  if (nrow(edges_dt) == 0L)
    return(data.frame(n_modules          = NA_integer_,
                      grey_rate          = NA_real_,
                      median_module_size = NA_real_))

  g  <- igraph::graph_from_data_frame(
    data.frame(gene_id_A = edges_dt$gene_id_A,
               gene_id_B = edges_dt$gene_id_B,
               weight    = edges_dt$mean_abs_r,
               stringsAsFactors = FALSE),
    directed = FALSE
  )
  set.seed(seed)
  cl    <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  sizes <- as.integer(igraph::sizes(cl))

  data.frame(
    n_modules          = as.integer(length(sizes)),
    grey_rate          = mean(sizes == 1L),
    median_module_size = stats::median(sizes)
  )
}
