library(testthat)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# AT-ID gene names
.gene_ids <- function(n) paste0("AT1G", formatC(seq_len(n), width = 5L, flag = "0"))

# Build a minimal InputBundle without going through Seurat.
# - n_genes: number of genes
# - n_cells: total cells; split evenly between strata
# - strata: character vector of stratum level names
# - groups_per_stratum: number of pseudobulk groups per stratum
# - seed: for reproducible counts
.make_bundle <- function(n_genes            = 50L,
                         n_cells_per_stratum = 10L,
                         strata             = c("Mock", "DC3000"),
                         groups_per_stratum = 5L,
                         seed               = 1L) {
  set.seed(seed)
  n_strata <- length(strata)
  n_cells  <- n_cells_per_stratum * n_strata
  gene_ids <- .gene_ids(n_genes)

  counts <- matrix(
    abs(rnorm(n_genes * n_cells, mean = 2, sd = 1)),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(gene_ids, paste0("cell_", seq_len(n_cells)))
  )

  stratum_vals <- rep(strata, each = n_cells_per_stratum)
  group_vals   <- paste0(
    "G",
    unlist(lapply(seq_len(n_strata), function(s) {
      offset <- (s - 1L) * groups_per_stratum
      rep(offset + seq_len(groups_per_stratum),
          length.out = n_cells_per_stratum)
    }))
  )

  cell_meta <- data.frame(
    cell_id   = paste0("cell_", seq_len(n_cells)),
    condition = stratum_vals,
    group_var = group_vals,
    stringsAsFactors = FALSE
  )

  gene_meta <- data.frame(
    gene_id     = gene_ids,
    gene_symbol = paste0("sym_", gene_ids),
    stringsAsFactors = FALSE
  )

  list(
    counts       = counts,
    cell_meta    = cell_meta,
    gene_meta    = gene_meta,
    stratum_spec = list(variable = "condition", levels = strata),
    dataset_id   = "test_dataset"
  )
}

# Predicate: does an object look like a valid NetworkResult?
.is_network_result <- function(x) {
  is.list(x) &&
    all(c("edge_table", "gene_ids", "stratum_id", "mode",
          "params", "timestamp") %in% names(x)) &&
    is.data.frame(x$edge_table) &&
    is.character(x$gene_ids) &&
    is.character(x$stratum_id) &&
    is.character(x$mode) &&
    is.list(x$params) &&
    inherits(x$timestamp, "POSIXct")
}

# ---------------------------------------------------------------------------
# estimate_pseudobulk tests
# ---------------------------------------------------------------------------

test_that("pseudobulk: returns a named list with names = stratum levels", {
  bundle <- .make_bundle(groups_per_stratum = 5L)
  res    <- estimate_pseudobulk(bundle, min_samples = 5L)

  expect_type(res, "list")
  expect_named(res, c("Mock", "DC3000"), ignore.order = TRUE)
})

test_that("pseudobulk: each element is a valid NetworkResult (all 6 slots, correct types)", {
  bundle <- .make_bundle(groups_per_stratum = 5L)
  res    <- estimate_pseudobulk(bundle, min_samples = 5L)

  for (nr in res) {
    expect_true(.is_network_result(nr),
                label = paste("NetworkResult valid for", nr$stratum_id))
  }
})

test_that("pseudobulk: edge_table has columns gene_id_A, gene_id_B, weight (numeric)", {
  bundle <- .make_bundle(groups_per_stratum = 5L)
  res    <- estimate_pseudobulk(bundle, min_samples = 5L)

  for (nr in res) {
    et <- nr$edge_table
    expect_true(all(c("gene_id_A", "gene_id_B", "weight") %in% names(et)),
                label = paste("edge_table columns for", nr$stratum_id))
    expect_type(et$weight, "double")
  }
})

test_that("pseudobulk: no self-pairs in edge_table", {
  bundle <- .make_bundle(groups_per_stratum = 5L)
  res    <- estimate_pseudobulk(bundle, min_samples = 5L)

  for (nr in res) {
    if (nrow(nr$edge_table) > 0L) {
      expect_true(all(nr$edge_table$gene_id_A != nr$edge_table$gene_id_B),
                  label = paste("no self-pairs in", nr$stratum_id))
    }
  }
})

test_that("pseudobulk: stratum with < min_samples pseudobulk groups is skipped with warning", {
  # Mock: 5 groups (passes min_samples = 5)
  # DC3000: 2 groups (skipped)
  bundle <- .make_bundle(n_cells_per_stratum = 10L, strata = c("Mock", "DC3000"))
  # Reassign DC3000 cells to only 2 groups
  is_dc <- bundle$cell_meta$condition == "DC3000"
  bundle$cell_meta$group_var[is_dc] <- rep(c("GX", "GY"), length.out = sum(is_dc))

  expect_warning(
    res <- estimate_pseudobulk(bundle, min_samples = 5L),
    regexp = "DC3000"
  )
  expect_true("Mock" %in% names(res))
  expect_false("DC3000" %in% names(res))
})

test_that("pseudobulk: params records group_var, n_cells_per_stratum, n_pseudobulk_per_stratum", {
  bundle <- .make_bundle(n_cells_per_stratum = 10L, groups_per_stratum = 5L)
  res    <- estimate_pseudobulk(bundle, min_samples = 5L)

  for (nr in res) {
    p <- nr$params
    expect_true("group_var" %in% names(p),
                label = paste("group_var in params of", nr$stratum_id))
    expect_equal(p$group_var, "group_var")

    expect_true("n_cells_per_stratum" %in% names(p))
    expect_type(p$n_cells_per_stratum, "integer")
    expect_named(p$n_cells_per_stratum, c("Mock", "DC3000"), ignore.order = TRUE)

    expect_true("n_pseudobulk_per_stratum" %in% names(p))
    expect_type(p$n_pseudobulk_per_stratum, "integer")
    expect_named(p$n_pseudobulk_per_stratum, c("Mock", "DC3000"), ignore.order = TRUE)
    expect_equal(p$n_pseudobulk_per_stratum["Mock"], c(Mock = 5L))
  }
})

test_that("pseudobulk: mode slot is 'pseudobulk'", {
  bundle <- .make_bundle(groups_per_stratum = 5L)
  res    <- estimate_pseudobulk(bundle, min_samples = 5L)
  for (nr in res) expect_equal(nr$mode, "pseudobulk")
})

# ---------------------------------------------------------------------------
# estimate_singlecellggm tests
# ---------------------------------------------------------------------------

test_that("singlecellggm: returns a named list with names = stratum levels", {
  skip_if_not_installed("corpcor")
  bundle <- .make_bundle(n_genes = 50L, n_cells_per_stratum = 10L)
  res    <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                   pcor_cutoff = 0.01, min_cells = 1L)

  expect_type(res, "list")
  expect_named(res, c("Mock", "DC3000"), ignore.order = TRUE)
})

test_that("singlecellggm: each element is a valid NetworkResult (all 6 slots)", {
  skip_if_not_installed("corpcor")
  bundle <- .make_bundle(n_genes = 50L, n_cells_per_stratum = 10L)
  res    <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                   pcor_cutoff = 0.01, min_cells = 1L)

  for (nr in res) {
    expect_true(.is_network_result(nr),
                label = paste("NetworkResult valid for", nr$stratum_id))
  }
})

test_that("singlecellggm: edge_table has columns gene_id_A, gene_id_B, weight", {
  skip_if_not_installed("corpcor")
  bundle <- .make_bundle(n_genes = 50L, n_cells_per_stratum = 10L)
  res    <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                   pcor_cutoff = 0.01, min_cells = 1L)

  for (nr in res) {
    et <- nr$edge_table
    expect_true(all(c("gene_id_A", "gene_id_B", "weight") %in% names(et)),
                label = paste("edge_table columns for", nr$stratum_id))
    if (nrow(et) > 0L) expect_type(et$weight, "double")
  }
})

test_that("singlecellggm: all gene IDs in edge_table are in bundle$gene_meta$gene_id", {
  skip_if_not_installed("corpcor")
  bundle <- .make_bundle(n_genes = 50L, n_cells_per_stratum = 10L)
  res    <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                   pcor_cutoff = 0.01, min_cells = 1L)

  valid_ids <- bundle$gene_meta$gene_id
  for (nr in res) {
    et <- nr$edge_table
    if (nrow(et) > 0L) {
      expect_true(all(et$gene_id_A %in% valid_ids),
                  label = paste("gene_id_A in gene_meta for", nr$stratum_id))
      expect_true(all(et$gene_id_B %in% valid_ids),
                  label = paste("gene_id_B in gene_meta for", nr$stratum_id))
    }
  }
})

test_that("singlecellggm: mode == 'singlecellggm'", {
  skip_if_not_installed("corpcor")
  bundle <- .make_bundle(n_genes = 50L, n_cells_per_stratum = 10L)
  res    <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                   pcor_cutoff = 0.01, min_cells = 1L)
  for (nr in res) expect_equal(nr$mode, "singlecellggm")
})

test_that("singlecellggm: params$aggregation == 'min_abs_pcor_across_iterations'", {
  skip_if_not_installed("corpcor")
  bundle <- .make_bundle(n_genes = 50L, n_cells_per_stratum = 10L)
  res    <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                   pcor_cutoff = 0.01, min_cells = 1L)
  for (nr in res) {
    expect_equal(nr$params$aggregation, "min_abs_pcor_across_iterations")
  }
})

test_that("singlecellggm: same seed produces identical edge_table on rerun", {
  skip_if_not_installed("corpcor")
  bundle <- .make_bundle(n_genes = 50L, n_cells_per_stratum = 10L)

  res1 <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                  pcor_cutoff = 0.01, min_cells = 1L, seed = 42L)
  res2 <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 20L,
                                  pcor_cutoff = 0.01, min_cells = 1L, seed = 42L)

  for (lvl in names(res1)) {
    et1 <- res1[[lvl]]$edge_table[order(res1[[lvl]]$edge_table$gene_id_A,
                                        res1[[lvl]]$edge_table$gene_id_B), ]
    et2 <- res2[[lvl]]$edge_table[order(res2[[lvl]]$edge_table$gene_id_A,
                                        res2[[lvl]]$edge_table$gene_id_B), ]
    expect_equal(et1, et2,
                 label = paste("reproducible edge_table for", lvl))
  }
})
