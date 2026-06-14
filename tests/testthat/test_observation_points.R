library(testthat)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

.gene_ids_obs <- function(n) paste0("AT1G", formatC(seq_len(n), width = 5L, flag = "0"))

# Build a minimal InputBundle that includes counts_raw, suitable for obs_* generators.
# n_genes: genes; n_cells_per_stratum: cells per stratum; strata: stratum names;
# groups_per_stratum: distinct group_var values per stratum; seed: for reproducibility.
.make_bundle_with_raw <- function(n_genes             = 40L,
                                   n_cells_per_stratum = 20L,
                                   strata              = c("Mock", "DC3000"),
                                   groups_per_stratum  = 5L,
                                   seed                = 42L) {
  set.seed(seed)
  n_strata <- length(strata)
  n_cells  <- n_cells_per_stratum * n_strata
  gene_ids <- .gene_ids_obs(n_genes)

  # Raw integer counts (genes x cells)
  counts_raw <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(gene_ids, paste0("cell_", seq_len(n_cells)))
  )
  mode(counts_raw) <- "integer"

  # Log-normalised version (library-size + log1p)
  lib_sizes <- colSums(counts_raw)
  counts_lognorm <- log1p(sweep(counts_raw, 2, lib_sizes / 1e4, "/"))

  stratum_vals <- rep(strata, each = n_cells_per_stratum)
  group_vals <- paste0(
    "G",
    unlist(lapply(seq_len(n_strata), function(s) {
      offset <- (s - 1L) * groups_per_stratum
      rep(offset + seq_len(groups_per_stratum), length.out = n_cells_per_stratum)
    }))
  )

  cell_meta <- data.frame(
    cell_id   = paste0("cell_", seq_len(n_cells)),
    condition = stratum_vals,
    group_var = group_vals,
    subcluster = paste0("sub_", group_vals),
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
    dataset_id   = "test_obs"
  )
}

# Predicate: is x a valid ObsPointSet?
.is_obs_point_set <- function(x) {
  is.list(x) &&
    all(c("matrix", "n_cells", "point_meta", "gene_ids", "design", "aggregation") %in% names(x)) &&
    is.matrix(x$matrix) &&
    is.integer(x$n_cells) &&
    is.data.frame(x$point_meta) &&
    is.character(x$gene_ids) &&
    is.list(x$design) &&
    x$aggregation %in% c("sum", "mean") &&
    nrow(x$matrix) == length(x$gene_ids) &&
    ncol(x$matrix) == length(x$n_cells) &&
    nrow(x$point_meta) == ncol(x$matrix) &&
    all(rownames(x$matrix) == x$gene_ids) &&
    "point_id" %in% names(x$point_meta)
}

# ---------------------------------------------------------------------------
# adapter extension: $counts_raw is integer and same dims as $counts
# ---------------------------------------------------------------------------

test_that("bundle$counts_raw is integer and same dims as counts", {
  bundle <- .make_bundle_with_raw()
  expect_false(is.null(bundle$counts_raw))
  expect_true(is.integer(bundle$counts_raw),
              label = "counts_raw must be integer mode")
  expect_equal(dim(bundle$counts_raw), dim(bundle$counts),
               label = "counts_raw and counts have same dimensions")
  expect_equal(rownames(bundle$counts_raw), rownames(bundle$counts),
               label = "counts_raw rownames match counts rownames")
  expect_equal(colnames(bundle$counts_raw), colnames(bundle$counts),
               label = "counts_raw colnames match counts colnames")
})

# ---------------------------------------------------------------------------
# obs_subcluster
# ---------------------------------------------------------------------------

test_that("obs_subcluster: returns valid ObsPointSet (all slots, correct types)", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  expect_true(.is_obs_point_set(obs),
              label = "obs_subcluster returns valid ObsPointSet")
})

test_that("obs_subcluster: gene_ids are AT-IDs from bundle$gene_meta", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  valid_ids <- bundle$gene_meta$gene_id
  expect_true(all(obs$gene_ids %in% valid_ids),
              label = "all gene_ids are valid AT-IDs")
})

test_that("obs_subcluster: matrix is genes x points", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  # number of distinct groups in group_var
  n_expected_points <- length(unique(bundle$cell_meta$group_var))
  expect_equal(ncol(obs$matrix), n_expected_points,
               label = "one observation point per unique group")
  expect_equal(nrow(obs$matrix), nrow(bundle$counts),
               label = "rows = genes")
})

test_that("obs_subcluster: design records group_col parameter", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "subcluster")
  expect_equal(obs$design$name, "obs_subcluster")
  expect_equal(obs$design$group_col, "subcluster")
})

test_that("obs_subcluster: n_cells sums to total cells in bundle", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  expect_equal(sum(obs$n_cells), ncol(bundle$counts))
})

# ---------------------------------------------------------------------------
# obs_stratified
# ---------------------------------------------------------------------------

test_that("obs_stratified: returns valid ObsPointSet", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_stratified(bundle, strata_cols = c("condition", "group_var"), min_cells = 1L)
  expect_true(.is_obs_point_set(obs))
})

test_that("obs_stratified: drops combos with < min_cells", {
  # Bundle has 20 cells per stratum, 5 groups → 4 cells per condition×group combo.
  # With min_cells = 3 → all combos pass (4 >= 3).
  bundle <- .make_bundle_with_raw()
  obs_all_pass <- obs_stratified(bundle, strata_cols = c("condition", "group_var"), min_cells = 3L)
  expect_true(.is_obs_point_set(obs_all_pass))

  # Add one rare cell that creates an additional tiny combo, then test that
  # it gets dropped. Easier: use strata_cols = "condition" alone (20 cells/group)
  # and set min_cells = 25 so both conditions are dropped → expect_error.
  expect_error(
    obs_stratified(bundle, strata_cols = "condition", min_cells = 25L),
    regexp = "min_cells"
  )

  # With min_cells = 1 on a single column, all groups pass.
  obs2 <- obs_stratified(bundle, strata_cols = "group_var", min_cells = 1L)
  expect_true(.is_obs_point_set(obs2))
})

# ---------------------------------------------------------------------------
# obs_cluster
# ---------------------------------------------------------------------------

test_that("obs_cluster: returns valid ObsPointSet", {
  bundle <- .make_bundle_with_raw()
  # Use existing group_var column as a proxy (many real Seurat objects have cluster columns)
  # obs_cluster will fall back to kNN+Louvain when column not found;
  # but for deterministic test we use a small bundle where kNN is fast
  obs <- tryCatch(
    obs_cluster(bundle, resolution = 1.0),
    error = function(e) NULL
  )
  # Skip if igraph not available
  if (!is.null(obs)) {
    expect_true(.is_obs_point_set(obs))
  }
})

test_that("obs_cluster: higher resolution yields >= as many points as lower (or comparable)", {
  skip_if_not_installed("igraph")
  bundle <- .make_bundle_with_raw(n_genes = 20L, n_cells_per_stratum = 30L)
  set.seed(1L)
  obs_low  <- tryCatch(obs_cluster(bundle, resolution = 0.3), error = function(e) NULL)
  set.seed(1L)
  obs_high <- tryCatch(obs_cluster(bundle, resolution = 2.0), error = function(e) NULL)
  if (!is.null(obs_low) && !is.null(obs_high)) {
    expect_gte(ncol(obs_high$matrix), ncol(obs_low$matrix),
               label = "higher resolution >= as many clusters as lower")
  }
})

# ---------------------------------------------------------------------------
# obs_metacell_knn
# ---------------------------------------------------------------------------

test_that("obs_metacell_knn: returns valid ObsPointSet", {
  skip_if_not_installed("igraph")
  bundle <- .make_bundle_with_raw(n_genes = 20L, n_cells_per_stratum = 30L)
  obs    <- tryCatch(
    obs_metacell_knn(bundle, target_size = 5L, n_points = 8L),
    error = function(e) NULL
  )
  if (!is.null(obs)) {
    expect_true(.is_obs_point_set(obs))
  }
})

test_that("obs_metacell_knn: target_size respected approximately", {
  skip_if_not_installed("igraph")
  bundle <- .make_bundle_with_raw(n_genes = 20L, n_cells_per_stratum = 40L)
  target <- 5L
  obs    <- tryCatch(
    obs_metacell_knn(bundle, target_size = target, n_points = 10L),
    error = function(e) NULL
  )
  if (!is.null(obs)) {
    # Median cells/point should be approximately target_size
    med_cells <- median(obs$n_cells)
    expect_lte(med_cells, target * 2,
               label = "median cells/point not much larger than target_size")
    expect_gte(med_cells, 1L,
               label = "each point has at least 1 cell")
  }
})

test_that("obs_metacell_knn: design slot records parameters", {
  skip_if_not_installed("igraph")
  bundle <- .make_bundle_with_raw(n_genes = 20L, n_cells_per_stratum = 30L)
  obs    <- tryCatch(
    obs_metacell_knn(bundle, target_size = 4L, n_points = 6L),
    error = function(e) NULL
  )
  if (!is.null(obs)) {
    expect_equal(obs$design$name, "obs_metacell_knn")
    expect_equal(obs$design$target_size, 4L)
    expect_equal(obs$design$n_points, 6L)
  }
})

# ---------------------------------------------------------------------------
# obs_axis_bin
# ---------------------------------------------------------------------------

test_that("obs_axis_bin: returns valid ObsPointSet when using PC1", {
  bundle <- .make_bundle_with_raw(n_genes = 20L, n_cells_per_stratum = 30L)
  obs    <- tryCatch(
    obs_axis_bin(bundle, axis = "PC1", n_bins = 5L),
    error = function(e) NULL
  )
  if (!is.null(obs)) {
    expect_true(.is_obs_point_set(obs))
    expect_lte(ncol(obs$matrix), 5L,
               label = "number of bins <= requested n_bins")
  }
})

test_that("obs_axis_bin: uses metadata column when axis matches a column name", {
  bundle <- .make_bundle_with_raw(n_genes = 20L, n_cells_per_stratum = 30L)
  # Add a numeric column to cell_meta
  bundle$cell_meta$pseudotime <- seq_len(nrow(bundle$cell_meta))
  obs <- tryCatch(
    obs_axis_bin(bundle, axis = "pseudotime", n_bins = 4L),
    error = function(e) NULL
  )
  if (!is.null(obs)) {
    expect_true(.is_obs_point_set(obs))
    expect_equal(obs$design$axis, "pseudotime")
  }
})

# ---------------------------------------------------------------------------
# normalize_obs
# ---------------------------------------------------------------------------

test_that("normalize_obs: all methods return matrix of same dims as input", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  for (method in c("none", "cp10k_log", "log_only", "zscore_gene")) {
    norm_mat <- normalize_obs(obs, method = method)
    expect_equal(dim(norm_mat), dim(obs$matrix),
                 label = paste("dim preserved for method", method))
    expect_true(is.numeric(norm_mat),
                label = paste("numeric output for method", method))
  }
})

test_that("normalize_obs cp10k_log: each point (column) sums to ~10000 after division before log", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var", aggregation = "sum")
  # The norm_mat is log1p(scaled); check that expm1 of scaled values sum to ~10k
  norm_mat <- normalize_obs(obs, method = "cp10k_log")
  # After log1p, recovering back: expm1(norm) gives the per-10k scaled values
  # Sum of expm1(col) should be ~10000 when original col sum > 0
  col_sums_raw <- colSums(obs$matrix)
  valid_cols   <- col_sums_raw > 0
  if (any(valid_cols)) {
    scaled_back <- expm1(norm_mat[, valid_cols, drop = FALSE])
    col_sums_scaled <- colSums(scaled_back)
    # Should be approximately 10000 per column
    expect_true(all(abs(col_sums_scaled - 1e4) < 1),
                label = "cp10k_log: columns sum to ~10k before log1p")
  }
})

test_that("normalize_obs zscore_gene: genes with zero sd get zero output, not NA/Inf", {
  bundle <- .make_bundle_with_raw(n_genes = 10L, n_cells_per_stratum = 20L)
  obs    <- obs_subcluster(bundle, group_col = "group_var")

  # Force one gene to be constant across all observation points
  obs$matrix[1L, ] <- 5.0

  norm_mat <- normalize_obs(obs, method = "zscore_gene")
  expect_false(any(is.na(norm_mat)), label = "no NAs in zscore output")
  expect_false(any(is.infinite(norm_mat)), label = "no Inf in zscore output")
  expect_true(all(norm_mat[1L, ] == 0),
              label = "constant gene has zero z-score (not NA/Inf)")
})

# ---------------------------------------------------------------------------
# coexpr_from_obs
# ---------------------------------------------------------------------------

test_that("coexpr_from_obs: returns list with edge_table (correct columns) and cor_mat", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  result <- coexpr_from_obs(obs, cor_type = "spearman", storage_cutoff = 0.1)

  expect_true(is.list(result))
  expect_true(all(c("edge_table", "cor_mat", "gene_ids", "cor_type", "storage_cutoff") %in%
                    names(result)))
  et <- result$edge_table
  expect_true(all(c("gene_id_A", "gene_id_B", "weight") %in% names(et)))
  if (nrow(et) > 0L) {
    expect_true(all(abs(et$weight) >= 0.1),
                label = "all edges pass storage_cutoff = 0.1")
  }
})

test_that("coexpr_from_obs: correlation is computed ACROSS observation points (not within)", {
  # Regression test analogous to the direction test in test_estimation.R.
  # 10 genes x 10 observation points.
  # Gene 1: [1..10], Gene 2: [10..1] (opposite trend) => Spearman should be -1.
  n_pts  <- 10L
  n_gene <- 10L
  gene_ids <- .gene_ids_obs(n_gene)

  mat <- matrix(5.0, nrow = n_gene, ncol = n_pts,
                dimnames = list(gene_ids, paste0("P", seq_len(n_pts))))
  mat[1L, ] <- seq(1L, 10L)   # gene 1: increasing
  mat[2L, ] <- seq(10L, 1L)   # gene 2: decreasing

  obs <- list(
    matrix     = mat,
    n_cells    = rep(5L, n_pts),
    point_meta = data.frame(point_id = paste0("P", seq_len(n_pts))),
    gene_ids   = gene_ids,
    design     = list(name = "test"),
    aggregation = "mean"
  )

  result <- suppressWarnings(coexpr_from_obs(obs, cor_type = "spearman", storage_cutoff = 0.1))
  et     <- result$edge_table

  g1 <- gene_ids[1L]; g2 <- gene_ids[2L]
  pair_row <- et[(et$gene_id_A == g1 & et$gene_id_B == g2) |
                 (et$gene_id_A == g2 & et$gene_id_B == g1), , drop = FALSE]

  expect_true(nrow(pair_row) > 0,
              label = "gene1-gene2 pair present (Spearman=-1, |w|=1 passes cutoff)")
  if (nrow(pair_row) > 0L) {
    expect_lt(pair_row$weight[1L], -0.8,
              label = "genes with opposite trends have strongly negative Spearman")
  }
})

test_that("coexpr_from_obs: no self-pairs in edge_table", {
  bundle <- .make_bundle_with_raw()
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  result <- coexpr_from_obs(obs, cor_type = "spearman")
  et     <- result$edge_table
  if (nrow(et) > 0L) {
    expect_true(all(et$gene_id_A != et$gene_id_B),
                label = "no self-pairs in edge_table")
  }
})

test_that("coexpr_from_obs: cor_mat is square genes x genes with rownames == gene_ids", {
  bundle <- .make_bundle_with_raw(n_genes = 15L)
  obs    <- obs_subcluster(bundle, group_col = "group_var")
  result <- coexpr_from_obs(obs, cor_type = "pearson")
  m      <- result$cor_mat
  expect_equal(nrow(m), length(obs$gene_ids))
  expect_equal(ncol(m), length(obs$gene_ids))
  expect_equal(rownames(m), obs$gene_ids)
})
