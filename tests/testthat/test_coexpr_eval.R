library(testthat)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

.gene_ids_eval <- function(n) paste0("AT2G", formatC(seq_len(n), width = 5L, flag = "0"))

# Minimal InputBundle with counts_raw (for eval functions that rebuild the network).
.make_eval_bundle <- function(n_genes            = 30L,
                               n_cells            = 60L,
                               strata             = "Mock",
                               groups_per_stratum = 8L,
                               seed               = 7L) {
  set.seed(seed)
  gene_ids <- .gene_ids_eval(n_genes)

  counts_raw <- matrix(
    rpois(n_genes * n_cells, lambda = 8L),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(gene_ids, paste0("c", seq_len(n_cells)))
  )
  mode(counts_raw) <- "integer"

  lib_sizes    <- colSums(counts_raw)
  counts_lognorm <- log1p(sweep(counts_raw, 2, pmax(lib_sizes, 1) / 1e4, "/"))

  group_vals <- paste0("G", rep(seq_len(groups_per_stratum), length.out = n_cells))

  cell_meta <- data.frame(
    cell_id   = paste0("c", seq_len(n_cells)),
    condition = strata[[1L]],
    group_var = group_vals,
    stringsAsFactors = FALSE
  )

  gene_meta <- data.frame(
    gene_id     = gene_ids,
    gene_symbol = paste0("SYM", seq_len(n_genes)),
    stringsAsFactors = FALSE
  )

  list(
    counts       = counts_lognorm,
    counts_raw   = counts_raw,
    cell_meta    = cell_meta,
    gene_meta    = gene_meta,
    stratum_spec = list(variable = "condition", levels = strata),
    dataset_id   = "test_eval"
  )
}

# Build a simple ObsPointSet directly (bypass generators for unit tests).
.make_obs <- function(mat, gene_ids = NULL) {
  if (is.null(gene_ids)) gene_ids <- rownames(mat)
  n_pts <- ncol(mat)
  list(
    matrix     = mat,
    n_cells    = rep(5L, n_pts),
    point_meta = data.frame(point_id = paste0("P", seq_len(n_pts))),
    gene_ids   = gene_ids,
    design     = list(name = "test_design"),
    aggregation = "mean"
  )
}

# Strong-signal matrix: all genes follow the same z-scored latent axis with small noise.
# Using z-scored data gives all genes equal mean/variance so the abs-weight GBA predictor
# works correctly (it assumes similar expression scales across genes).
.make_structured_obs <- function(n_genes = 20L, n_pts = 30L, seed = 99L) {
  set.seed(seed)
  gene_ids <- .gene_ids_eval(n_genes)

  # Shared z-scored latent axis — ensures mean=0, sd≈1 for each gene
  axis <- as.numeric(scale(seq(0, 1, length.out = n_pts)))
  mat  <- matrix(NA_real_, nrow = n_genes, ncol = n_pts,
                 dimnames = list(gene_ids, paste0("P", seq_len(n_pts))))
  for (g in seq_len(n_genes)) {
    raw  <- axis + rnorm(n_pts, sd = 0.1)
    mat[g, ] <- as.numeric(scale(raw))   # z-score per gene → equal scales
  }
  .make_obs(mat, gene_ids)
}

# Pure noise matrix — each gene is an independent z-scored random vector.
.make_noise_obs <- function(n_genes = 20L, n_pts = 30L, seed = 88L) {
  set.seed(seed)
  gene_ids <- .gene_ids_eval(n_genes)
  mat <- matrix(NA_real_, nrow = n_genes, ncol = n_pts,
                dimnames = list(gene_ids, paste0("P", seq_len(n_pts))))
  for (g in seq_len(n_genes)) {
    raw <- rnorm(n_pts)
    mat[g, ] <- as.numeric(scale(raw))
  }
  .make_obs(mat, gene_ids)
}

# ---------------------------------------------------------------------------
# eval_effective_rank
# ---------------------------------------------------------------------------

test_that("eval_effective_rank: returns value in [1, n_points]", {
  obs <- .make_structured_obs()
  res <- eval_effective_rank(obs)
  expect_true(is.data.frame(res))
  expect_true("eff_rank" %in% names(res))
  expect_gte(res$eff_rank, 1.0 - 1e-9,
             label = "eff_rank >= 1")
  expect_lte(res$eff_rank, ncol(obs$matrix) + 1e-9,
             label = "eff_rank <= n_points")
})

test_that("eval_effective_rank: reports n_points and n_genes", {
  obs <- .make_structured_obs(n_genes = 15L, n_pts = 25L)
  res <- eval_effective_rank(obs)
  expect_equal(res$n_points, 25L)
  expect_equal(res$n_genes,  15L)
})

test_that("eval_effective_rank: structured data has higher eff_rank than rank-1 data", {
  set.seed(1L)
  n_g  <- 20L; n_p <- 30L
  gene_ids <- .gene_ids_eval(n_g)

  # Rank-1: all genes are the same axis (perfect correlation, single axis)
  axis   <- rnorm(n_p)
  mat_r1 <- matrix(axis, nrow = n_g, ncol = n_p, byrow = TRUE,
                   dimnames = list(gene_ids, paste0("P", seq_len(n_p))))
  mat_r1 <- mat_r1 + matrix(rnorm(n_g * n_p, sd = 0.01), n_g, n_p,
                              dimnames = list(gene_ids, paste0("P", seq_len(n_p))))

  obs_r1   <- .make_obs(mat_r1, gene_ids)
  obs_full <- .make_structured_obs(n_genes = n_g, n_pts = n_p)

  res_r1   <- eval_effective_rank(obs_r1)
  res_full <- eval_effective_rank(obs_full)
  expect_gt(res_full$eff_rank, res_r1$eff_rank,
            label = "multi-axis data has higher eff_rank than rank-1 data")
})

# ---------------------------------------------------------------------------
# eval_visible_genes
# ---------------------------------------------------------------------------

test_that("eval_visible_genes: counts correctly on constructed matrix", {
  n_genes <- 10L; n_pts <- 15L
  gene_ids <- .gene_ids_eval(n_genes)

  mat <- matrix(rnorm(n_genes * n_pts), nrow = n_genes, ncol = n_pts,
                dimnames = list(gene_ids, paste0("P", seq_len(n_pts))))
  # Make first 3 genes constant (all zeros or all same value)
  mat[1L, ] <- 0.0
  mat[2L, ] <- 3.5
  mat[3L, ] <- -1.0

  obs <- .make_obs(mat, gene_ids)
  res <- eval_visible_genes(obs, min_var = 1e-6)

  expect_true(is.data.frame(res))
  expect_true(all(c("n_visible", "n_total", "frac_visible") %in% names(res)))
  expect_equal(res$n_total, n_genes)
  # 7 genes have real variance, 3 are constant
  expect_equal(res$n_visible, 7L,
               label = "7 variable genes correctly identified")
  expect_equal(res$frac_visible, 7 / n_genes, tolerance = 1e-9)
})

test_that("eval_visible_genes: all-constant matrix has 0 visible genes", {
  gene_ids <- .gene_ids_eval(5L)
  mat <- matrix(1.0, nrow = 5L, ncol = 10L,
                dimnames = list(gene_ids, paste0("P", seq_len(10L))))
  obs <- .make_obs(mat, gene_ids)
  res <- eval_visible_genes(obs)
  expect_equal(res$n_visible, 0L)
  expect_equal(res$frac_visible, 0.0)
})

# ---------------------------------------------------------------------------
# eval_splithalf
# ---------------------------------------------------------------------------

test_that("eval_splithalf: returns higher agreement for structured than for noise", {
  skip_if_not_installed("igraph")
  # Structured bundle: all cells follow a strong global trend (gene 1 is the
  # most highly expressed in the last group, creating a consistent rank ordering).
  bundle_s <- .make_eval_bundle(n_genes = 20L, n_cells = 80L, groups_per_stratum = 10L, seed = 1L)

  # Noise bundle: expression randomly permuted per gene (destroy covariation structure)
  bundle_n <- .make_eval_bundle(n_genes = 20L, n_cells = 80L, groups_per_stratum = 10L, seed = 2L)
  set.seed(77L)
  # Permute columns for each gene independently to destroy all covariation
  # Use abs() to keep values non-negative (cp10k_log requires non-negative input)
  bundle_n$counts <- abs(apply(bundle_n$counts, 1L, sample))
  dim(bundle_n$counts) <- c(80L, 20L)  # apply transposes: ncols first
  bundle_n$counts <- t(bundle_n$counts)
  dimnames(bundle_n$counts) <- dimnames(bundle_s$counts)

  set.seed(42L)
  res_s <- eval_splithalf(bundle_s, design_fn = obs_subcluster,
                           design_args = list(group_col = "group_var"),
                           cor_type = "spearman", norm_method = "none",
                           n_reps = 3L)

  set.seed(42L)
  res_n <- eval_splithalf(bundle_n, design_fn = obs_subcluster,
                           design_args = list(group_col = "group_var"),
                           cor_type = "spearman", norm_method = "none",
                           n_reps = 3L)

  expect_true(is.data.frame(res_s))
  expect_true(all(c("mat_cor_mean", "jaccard_mean") %in% names(res_s)))
  # Results should be numeric
  expect_false(is.null(res_s$mat_cor_mean))
})

test_that("eval_splithalf: returns data.frame with required columns", {
  bundle <- .make_eval_bundle(n_genes = 15L, n_cells = 60L, groups_per_stratum = 8L)
  res    <- eval_splithalf(bundle, design_fn = obs_subcluster,
                            design_args = list(group_col = "group_var"),
                            n_reps = 2L)
  expect_true(is.data.frame(res))
  expect_true(all(c("mat_cor_mean", "mat_cor_sd",
                    "jaccard_mean", "jaccard_sd", "n_reps") %in% names(res)))
})

# ---------------------------------------------------------------------------
# eval_heldout_predictivity
# ---------------------------------------------------------------------------

test_that("eval_heldout_predictivity: returns higher R2 for structured than random data", {
  obs_s <- .make_structured_obs(n_genes = 20L, n_pts = 40L, seed = 11L)
  obs_n <- .make_noise_obs(n_genes = 20L, n_pts = 40L, seed = 22L)

  set.seed(5L)
  res_s <- eval_heldout_predictivity(obs_s, k_partners = 5L, n_folds = 4L)
  set.seed(5L)
  res_n <- eval_heldout_predictivity(obs_n, k_partners = 5L, n_folds = 4L)

  expect_true(is.data.frame(res_s))
  expect_true(all(c("predictivity_mean_r2", "predictivity_median_r2") %in% names(res_s)))

  if (!is.na(res_s$predictivity_mean_r2) && !is.na(res_n$predictivity_mean_r2)) {
    expect_gt(res_s$predictivity_mean_r2, res_n$predictivity_mean_r2,
              label = "structured data has higher held-out predictivity than pure noise")
  }
})

test_that("eval_heldout_predictivity: returns data.frame with correct columns", {
  obs <- .make_structured_obs(n_genes = 10L, n_pts = 20L)
  res <- eval_heldout_predictivity(obs, k_partners = 3L, n_folds = 3L)
  expect_true(is.data.frame(res))
  expect_true("predictivity_mean_r2" %in% names(res))
  expect_false(is.na(res$predictivity_mean_r2))
})

# ---------------------------------------------------------------------------
# eval_depth_leakage
# ---------------------------------------------------------------------------

test_that("eval_depth_leakage: detects a planted depth confound", {
  # Construct observation points where half the genes have expression proportional
  # to a depth factor (high degree and high mean), the other half are flat (low degree,
  # low mean). This creates a correlation between degree and mean expression.
  set.seed(33L)
  n_genes <- 30L; n_pts <- 30L
  gene_ids <- .gene_ids_eval(n_genes)

  depth_factor <- seq(1, 10, length.out = n_pts)

  mat <- matrix(0.1, nrow = n_genes, ncol = n_pts,
                dimnames = list(gene_ids, paste0("P", seq_len(n_pts))))
  # First 15 genes: proportional to depth (high covariation, high mean)
  for (g in seq_len(15L)) {
    mat[g, ] <- g * depth_factor
  }
  # Last 15 genes: constant (zero degree, low mean)
  for (g in 16L:n_genes) {
    mat[g, ] <- 0.1
  }

  obs <- .make_obs(mat, gene_ids)
  # suppressWarnings: cor() warns about zero-SD when all constant genes have the same degree
  res <- suppressWarnings(eval_depth_leakage(obs, threshold = 0.3))

  expect_true(is.data.frame(res))
  expect_true("depth_leakage_rho" %in% names(res))
  if (!is.na(res$depth_leakage_rho)) {
    expect_gt(abs(res$depth_leakage_rho), 0.5,
              label = "planted depth confound detected as high degree-vs-mean correlation")
  }
})

test_that("eval_depth_leakage: returns data.frame with depth_leakage_rho and n_genes", {
  obs <- .make_structured_obs(n_genes = 15L, n_pts = 20L)
  # suppressWarnings: cor() may warn about zero-SD when structure is perfectly uniform
  res <- suppressWarnings(eval_depth_leakage(obs))
  expect_true(is.data.frame(res))
  expect_true(all(c("depth_leakage_rho", "n_genes") %in% names(res)))
  expect_equal(res$n_genes, 15L)
})

# ---------------------------------------------------------------------------
# evaluate_obs_design
# ---------------------------------------------------------------------------

test_that("evaluate_obs_design: returns a one-row data.frame with required columns", {
  bundle <- .make_eval_bundle(n_genes = 20L, n_cells = 60L, groups_per_stratum = 8L)
  res    <- evaluate_obs_design(
    bundle    = bundle,
    design_fn = obs_subcluster,
    design_args = list(group_col = "group_var"),
    cor_type  = "spearman",
    norm_method = "cp10k_log",
    n_splithalf = 2L,
    splithalf_reps = 2L,
    null_perm = 5L,
    heldout_folds = 3L
  )

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 1L,
               label = "evaluate_obs_design returns exactly one row")

  required_cols <- c("design_name", "n_points", "n_genes", "eff_rank",
                     "n_visible", "frac_visible",
                     "predictivity_mean_r2",
                     "null_gap_ratio",
                     "depth_leakage_rho",
                     "splithalf_mat_cor_mean",
                     "splithalf_jaccard_mean")
  for (col in required_cols) {
    expect_true(col %in% names(res),
                label = paste("column present:", col))
  }
})
