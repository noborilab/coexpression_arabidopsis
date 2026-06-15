#' @title Stage 3 Threshold Evaluation Adapters
#'
#' @description
#' Thin adapters that allow the prior-free metrics from `coexpr_eval.R` to
#' evaluate THRESHOLDED NETWORKS rather than ObsPointSet / InputBundle objects.
#'
#' **Stage 3 context:** The 4 pathogen conditions (Mock / DC3000 / AvrRpt2 /
#' AvrRpm1) play the role of "observation points." Each adapter builds a
#' **genes × 4-conditions fingerprint matrix** where entry [g, c] is the mean
#' signed Spearman r between gene g and all its retained network neighbours in
#' condition c. This matrix is passed to the corresponding `coexpr_eval.R`
#' function unchanged; only the *input construction* differs from the
#' ObsPointSet path.
#'
#' **Metric interpretations in Stage 3:**
#' - `eval_splithalf`  → Jaccard stability when conditions split 2+2 (3 unique
#'   complementary splits); each half's retained edges are those above the
#'   threshold / in the top-k for the half-specific mean |r|.
#' - `eval_effective_rank` → participation ratio of singular values of the
#'   fingerprint matrix; higher = network fingerprints vary along more
#'   independent axes across conditions (≤ 4 by construction with 4 conditions).
#' - `eval_null_gap`   → fraction of gene-pair fingerprint correlations above
#'   0.3 in real vs. row-permuted fingerprint matrix (gene subsample ≤ 2000).
#' - `eval_heldout_predictivity` → leave-one-condition-out GBA R² on the
#'   fingerprint matrix (4-fold CV, 500-gene subsample guard preserved).
#' - `eval_visible_genes` → non-isolated genes with ≥ 1 retained edge.
#'
#' **What is NOT adapted (and why):**
#' - `eval_downsample_depth`  — requires raw count matrix; N/A.
#' - `eval_downsample_cells`  — requires raw cell matrix; N/A.
#' - `eval_depth_leakage`     — depth confounding is pre-adjusted; not a
#'   threshold design dimension.
#'
#' @name stage3_threshold_eval
NULL

# ---------------------------------------------------------------------------
# Internal helper: build genes × conditions fingerprint matrix
# ---------------------------------------------------------------------------

#' Build a genes × conditions fingerprint matrix
#'
#' For each gene g and condition c the fingerprint value is the mean signed
#' Spearman r between g and all its retained network neighbours in condition c.
#' Computed from `edges_dt` which must contain per-condition r columns.
#'
#' @param edges_dt data.table with columns `gene_id_A`, `gene_id_B`, and one
#'   numeric column per condition in `cor_cols`.
#' @param cor_cols Character vector naming the per-condition r columns.
#' @return List with `$matrix` (genes × conditions) and `$gene_ids`.
#' @keywords internal
.stage3_fingerprint <- function(edges_dt,
                                 cor_cols = c("r_Mock", "r_DC3000",
                                              "r_AvrRpt2", "r_AvrRpm1")) {
  if (nrow(edges_dt) == 0L) {
    m <- matrix(numeric(0), nrow = 0L, ncol = length(cor_cols),
                dimnames = list(NULL, cor_cols))
    return(list(matrix = m, gene_ids = character(0L),
                design = list(name = "stage3_network")))
  }

  genes_in_net <- sort(unique(c(edges_dt$gene_id_A, edges_dt$gene_id_B)))
  n_genes <- length(genes_in_net)
  n_conds <- length(cor_cols)

  mat <- matrix(0, nrow = n_genes, ncol = n_conds,
                dimnames = list(genes_in_net, cor_cols))
  cnt <- matrix(0L, nrow = n_genes, ncol = n_conds)

  gene_idx <- setNames(seq_len(n_genes), genes_in_net)
  idx_A    <- gene_idx[edges_dt$gene_id_A]
  idx_B    <- gene_idx[edges_dt$gene_id_B]

  for (ci in seq_len(n_conds)) {
    r_vals <- edges_dt[[cor_cols[ci]]]
    # Both directions: A's neighbour B and B's neighbour A
    for (grp in list(list(idx = idx_A, r = r_vals),
                     list(idx = idx_B, r = r_vals))) {
      ok      <- !is.na(grp$r)
      s_names <- tapply(grp$r[ok], grp$idx[ok], sum, na.rm = TRUE)
      c_tab   <- tabulate(grp$idx[ok], nbins = n_genes)
      row_ids <- as.integer(names(s_names))
      mat[row_ids, ci] <- mat[row_ids, ci] + s_names
      cnt[, ci]        <- cnt[, ci] + c_tab
    }
  }

  # Mean = sum / count; 0-count genes stay 0 (safe divisor via pmax)
  cnt_safe <- pmax(cnt, 1L)
  mat      <- mat / cnt_safe

  list(matrix = mat, gene_ids = genes_in_net,
       design = list(name = "stage3_network"))
}

# ---------------------------------------------------------------------------
# stage3_eval_visible_genes
# ---------------------------------------------------------------------------

#' Non-isolated gene count for a thresholded network
#'
#' @param edges_dt data.table with columns `gene_id_A` and `gene_id_B`.
#' @param n_total Integer; total gene universe size. Default `11010L`.
#' @return One-row data.frame: `n_visible`, `n_total`, `frac_visible`.
#' @export
stage3_eval_visible_genes <- function(edges_dt, n_total = 11010L) {
  if (nrow(edges_dt) == 0L)
    return(data.frame(n_visible = 0L, n_total = n_total, frac_visible = 0))

  n_visible <- length(unique(c(edges_dt$gene_id_A, edges_dt$gene_id_B)))
  data.frame(n_visible    = as.integer(n_visible),
             n_total      = as.integer(n_total),
             frac_visible = n_visible / n_total)
}

# ---------------------------------------------------------------------------
# stage3_eval_effective_rank
# ---------------------------------------------------------------------------

#' Effective rank of the genes × conditions fingerprint matrix
#'
#' Calls `eval_effective_rank` on the fingerprint. Maximum possible value is
#' `length(cor_cols)` (= 4 for four conditions).
#'
#' @inheritParams .stage3_fingerprint
#' @return One-row data.frame: `eff_rank`, `n_points`, `n_genes`.
#' @export
stage3_eval_effective_rank <- function(edges_dt,
                                        cor_cols = c("r_Mock", "r_DC3000",
                                                     "r_AvrRpt2", "r_AvrRpm1")) {
  obs <- .stage3_fingerprint(edges_dt, cor_cols)
  eval_effective_rank(obs)
}

# ---------------------------------------------------------------------------
# stage3_eval_null_gap
# ---------------------------------------------------------------------------

#' Null-gap ratio for the fingerprint matrix
#'
#' Calls `eval_null_gap` on the fingerprint matrix. To avoid OOM on dense
#' networks (genes × 4 fingerprint with > 2000 genes), subsamples to at most
#' `max_genes` genes before computing the n × n correlation matrix.
#'
#' The threshold (default 0.3) is applied to Spearman correlations between
#' gene fingerprints (not to the original edge |r|), testing whether gene-pair
#' fingerprints cluster together more than expected under a row-permuted null.
#'
#' @inheritParams .stage3_fingerprint
#' @param n_perm Integer; permutations for the null. Default `10L`.
#' @param threshold Fingerprint correlation threshold. Default `0.3`.
#' @param max_genes Integer; gene subsample ceiling. Default `2000L`.
#' @return One-row data.frame: `null_gap_ratio`, `real_frac`, `perm_frac_mean`.
#' @export
stage3_eval_null_gap <- function(edges_dt,
                                  cor_cols  = c("r_Mock", "r_DC3000",
                                                "r_AvrRpt2", "r_AvrRpm1"),
                                  n_perm    = 10L,
                                  threshold = 0.3,
                                  max_genes = 2000L) {
  obs <- .stage3_fingerprint(edges_dt, cor_cols)
  ng  <- nrow(obs$matrix)

  if (ng < 2L)
    return(data.frame(null_gap_ratio = NA_real_,
                      real_frac      = NA_real_,
                      perm_frac_mean = NA_real_))

  # Subsample to avoid O(n²) OOM with large networks
  if (ng > max_genes) {
    keep <- sort(sample(ng, max_genes))
    obs$matrix   <- obs$matrix[keep, , drop = FALSE]
    obs$gene_ids <- obs$gene_ids[keep]
    message("stage3_eval_null_gap: subsampled to ", max_genes,
            " genes (from ", ng, ") to bound O(n²) correlation matrix.")
  }

  eval_null_gap(obs, cor_type = "spearman", n_perm = n_perm,
                threshold = threshold)
}

# ---------------------------------------------------------------------------
# stage3_eval_heldout_predictivity
# ---------------------------------------------------------------------------

#' Held-out predictivity on the fingerprint matrix
#'
#' Calls `eval_heldout_predictivity` with leave-one-condition-out CV
#' (n_folds = 4). The 500-gene subsample guard in the original function is
#' preserved.
#'
#' @inheritParams .stage3_fingerprint
#' @param k_partners Integer; top-k partners for GBA prediction. Default `10L`.
#' @return One-row data.frame: `predictivity_mean_r2`, `predictivity_median_r2`.
#' @export
stage3_eval_heldout_predictivity <- function(edges_dt,
                                              cor_cols   = c("r_Mock",
                                                             "r_DC3000",
                                                             "r_AvrRpt2",
                                                             "r_AvrRpm1"),
                                              k_partners = 10L) {
  obs <- .stage3_fingerprint(edges_dt, cor_cols)
  if (nrow(obs$matrix) < 2L || ncol(obs$matrix) < 4L)
    return(data.frame(predictivity_mean_r2   = NA_real_,
                      predictivity_median_r2 = NA_real_))
  eval_heldout_predictivity(obs, k_partners = k_partners,
                            n_folds  = ncol(obs$matrix),  # 4 = leave-one-out
                            cor_type = "spearman")
}

# ---------------------------------------------------------------------------
# stage3_eval_splithalf
# ---------------------------------------------------------------------------

#' Condition-split stability of a thresholded network
#'
#' Splits the 4 conditions into 3 unique complementary 2+2 pairs and computes
#' Jaccard similarity of retained edges between the two halves. This tests
#' whether the same edges survive the threshold when computed on different
#' condition subsets.
#'
#' @param cand_dt data.table of candidate pairs with columns `gene_id_A`,
#'   `gene_id_B`, and per-condition r columns (`cor_cols`). Must include all
#'   pairs that plausibly survive the threshold in either half:
#'   - Lever A: pairs within 0.08 below the design threshold.
#'   - Lever B: top-200-per-gene candidates (gives ≥ top-100 buffer per half).
#' @param threshold_r Numeric scalar; global-|r| threshold (Lever A). Exactly
#'   one of `threshold_r` or `topk` must be set.
#' @param topk Integer; per-gene top-k (Lever B). Applied per-gene to the
#'   half-specific mean |r|; union of top-k for each endpoint.
#' @param cor_cols Character vector of per-condition r column names.
#' @return One-row data.frame: `splithalf_jaccard`, `splithalf_jaccard_sd`,
#'   `n_splits`.
#' @export
stage3_eval_splithalf <- function(cand_dt,
                                   threshold_r = NULL,
                                   topk        = NULL,
                                   cor_cols    = c("r_Mock", "r_DC3000",
                                                   "r_AvrRpt2", "r_AvrRpm1")) {
  if (is.null(threshold_r) == is.null(topk))
    stop("stage3_eval_splithalf: exactly one of threshold_r or topk must be set.")

  # 3 unique complementary 2+2 splits
  splits <- list(
    list(A = 1:2,      B = 3:4),
    list(A = c(1, 3),  B = c(2, 4)),
    list(A = c(1, 4),  B = c(2, 3))
  )

  jaccards <- vapply(splits, function(sp) {
    r_A <- rowMeans(abs(as.matrix(cand_dt[, cor_cols[sp$A], with = FALSE])),
                    na.rm = TRUE)
    r_B <- rowMeans(abs(as.matrix(cand_dt[, cor_cols[sp$B], with = FALSE])),
                    na.rm = TRUE)

    if (!is.null(threshold_r)) {
      in_A <- r_A >= threshold_r
      in_B <- r_B >= threshold_r
    } else {
      # Per-gene top-k; union rule: edge retained if in top-k for EITHER endpoint
      gA <- cand_dt$gene_id_A
      gB <- cand_dt$gene_id_B
      # ave: for each group gi, rank() applied to the sub-vector → lower rank = higher |r|
      rank_A_as_A <- ave(-r_A, gA, FUN = rank)
      rank_A_as_B <- ave(-r_A, gB, FUN = rank)
      rank_B_as_A <- ave(-r_B, gA, FUN = rank)
      rank_B_as_B <- ave(-r_B, gB, FUN = rank)

      in_A <- (rank_A_as_A <= topk) | (rank_A_as_B <= topk)
      in_B <- (rank_B_as_A <= topk) | (rank_B_as_B <= topk)
    }

    n_int <- sum(in_A & in_B, na.rm = TRUE)
    n_uni <- sum(in_A | in_B, na.rm = TRUE)
    if (n_uni == 0L) NA_real_ else n_int / n_uni
  }, numeric(1L))

  n_ok <- sum(!is.na(jaccards))
  data.frame(
    splithalf_jaccard    = mean(jaccards, na.rm = TRUE),
    splithalf_jaccard_sd = if (n_ok > 1L) stats::sd(jaccards, na.rm = TRUE)
                          else NA_real_,
    n_splits             = as.integer(n_ok)
  )
}

# ---------------------------------------------------------------------------
# stage3_eval_louvain
# ---------------------------------------------------------------------------

#' Fast Louvain module statistics for a thresholded network
#'
#' Runs one-pass Louvain community detection (igraph) on the undirected
#' network defined by `edges_dt`. Weights are absolute mean |r| values.
#' Intended as descriptive support, NOT a selection metric.
#'
#' @param edges_dt data.table with columns `gene_id_A`, `gene_id_B`, `abs_r`.
#' @param seed Integer RNG seed. Default `98L`.
#' @return One-row data.frame: `n_modules`, `grey_rate` (fraction of size-1
#'   modules), `median_module_size`.
#' @export
stage3_eval_louvain <- function(edges_dt, seed = 98L) {
  if (nrow(edges_dt) == 0L)
    return(data.frame(n_modules          = NA_integer_,
                      grey_rate          = NA_real_,
                      median_module_size = NA_real_))

  g  <- igraph::graph_from_data_frame(
    edges_dt[, .(gene_id_A, gene_id_B, weight = abs_r)],
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
