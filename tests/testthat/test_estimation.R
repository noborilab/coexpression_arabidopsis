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

# Build a bundle suitable for GGM tests (dense counts, many cells).
.make_ggm_bundle <- function(n_genes             = 40L,
                              n_cells_per_stratum = 200L,
                              strata              = c("Mock", "DC3000"),
                              seed                = 7L) {
  set.seed(seed)
  n_strata <- length(strata)
  n_cells  <- n_cells_per_stratum * n_strata
  gene_ids <- .gene_ids(n_genes)

  counts <- matrix(
    abs(rnorm(n_genes * n_cells, mean = 2, sd = 1)),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(gene_ids, paste0("cell_", seq_len(n_cells)))
  )

  cell_meta <- data.frame(
    cell_id   = paste0("cell_", seq_len(n_cells)),
    condition = rep(strata, each = n_cells_per_stratum),
    group_var = paste0("G", rep(seq_len(5L), length.out = n_cells)),
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
    dataset_id   = "test_ggm"
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
  bundle <- .make_ggm_bundle()
  res    <- estimate_singlecellggm(bundle, n_iter = 3L, subsample = 10L,
                                    pcor_cutoff = 0.0, keep_negative = TRUE)
  expect_type(res, "list")
  expect_named(res, c("Mock", "DC3000"), ignore.order = TRUE)
})

test_that("singlecellggm: each element is a valid NetworkResult (all 6 slots)", {
  bundle <- .make_ggm_bundle()
  res    <- estimate_singlecellggm(bundle, n_iter = 3L, subsample = 10L,
                                    pcor_cutoff = 0.0, keep_negative = TRUE)
  for (nr in res) {
    expect_true(.is_network_result(nr),
                label = paste("NetworkResult valid for", nr$stratum_id))
  }
})

test_that("singlecellggm: edge_table has expected columns including coex_cells and sampling_num", {
  bundle <- .make_ggm_bundle()
  res    <- estimate_singlecellggm(bundle, n_iter = 3L, subsample = 10L,
                                    pcor_cutoff = 0.0, keep_negative = TRUE)
  for (nr in res) {
    et <- nr$edge_table
    expect_true(
      all(c("gene_id_A", "gene_id_B", "weight", "coex_cells", "sampling_num") %in% names(et)),
      label = paste("edge_table columns for", nr$stratum_id)
    )
    if (nrow(et) > 0L) {
      expect_type(et$weight,       "double")
      expect_type(et$coex_cells,   "integer")
      expect_type(et$sampling_num, "integer")
    }
  }
})

test_that("singlecellggm: mode == 'singlecellggm'", {
  bundle <- .make_ggm_bundle()
  res    <- estimate_singlecellggm(bundle, n_iter = 3L, subsample = 10L)
  for (nr in res) expect_equal(nr$mode, "singlecellggm")
})

test_that("singlecellggm: params$aggregation == 'min_abs_pcor_across_iterations'", {
  bundle <- .make_ggm_bundle()
  res    <- estimate_singlecellggm(bundle, n_iter = 3L, subsample = 10L)
  for (nr in res) {
    expect_equal(nr$params$aggregation, "min_abs_pcor_across_iterations")
  }
})

test_that("singlecellggm: n_iter auto-resolves (NULL + p <= subsample → n_iter = 1)", {
  # p = 40, subsample = 2000 (default) → subsample >= p → n_iter = 1, subsample = p
  bundle <- .make_ggm_bundle(n_genes = 40L)
  res <- NULL
  expect_no_error({
    res <- estimate_singlecellggm(bundle, n_iter = NULL, subsample = 2000L)
  })
  for (nr in res) {
    expect_true(.is_network_result(nr))
    expect_equal(nr$params$n_iter, 1L)
  }
})

test_that("singlecellggm: n_iter auto-resolves formula when p > subsample", {
  # p = 40, subsample = 10 → n_iter = max(1, round(40*39/39980)) = 1
  bundle <- .make_ggm_bundle(n_genes = 40L)
  res    <- estimate_singlecellggm(bundle, n_iter = NULL, subsample = 10L)
  for (nr in res) {
    p_actual        <- nr$params$n_genes
    expected_n_iter <- max(1L, round(p_actual * (p_actual - 1L) / 39980))
    expect_equal(nr$params$n_iter, expected_n_iter)
  }
})

test_that("singlecellggm: same seed produces identical edge_table (reproducibility)", {
  bundle <- .make_ggm_bundle(n_genes = 40L, n_cells_per_stratum = 200L)

  res1 <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 15L,
                                  pcor_cutoff = 0.0, keep_negative = TRUE, seed = 42L)
  res2 <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 15L,
                                  pcor_cutoff = 0.0, keep_negative = TRUE, seed = 42L)

  for (lvl in names(res1)) {
    et1 <- res1[[lvl]]$edge_table
    et2 <- res2[[lvl]]$edge_table
    et1 <- et1[order(et1$gene_id_A, et1$gene_id_B), ]
    et2 <- et2[order(et2$gene_id_A, et2$gene_id_B), ]
    rownames(et1) <- rownames(et2) <- NULL
    expect_equal(et1, et2, label = paste("reproducible edge_table for", lvl))
  }
})

test_that("singlecellggm: coex filter excludes low-co-detection pairs", {
  # Gene 1: detected cells 1:15 (15 cells); Gene 2: detected cells 11:30 (20 cells)
  # Co-detection: cells 11:15 = 5 cells < coex_cutoff = 10 → pair excluded
  # Both genes still pass the per-gene filter (>= 10 detected cells each)
  set.seed(42L)
  n_genes <- 30L
  n_cells <- 200L
  gene_ids <- .gene_ids(n_genes)

  counts <- matrix(
    abs(rnorm(n_genes * n_cells, mean = 3, sd = 0.5)),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(gene_ids, paste0("c", seq_len(n_cells)))
  )
  counts[1L, 16:200]             <- 0   # gene 1: only cells 1:15
  counts[2L, c(1:10, 31:200)]   <- 0   # gene 2: only cells 11:30
  # co-detection of (gene1, gene2) = cells 11:15 = 5 cells

  bundle <- list(
    counts       = counts,
    cell_meta    = data.frame(cell_id   = paste0("c", seq_len(n_cells)),
                               condition = "Mock",
                               group_var = "G1",
                               stringsAsFactors = FALSE),
    gene_meta    = data.frame(gene_id     = gene_ids,
                               gene_symbol = gene_ids,
                               stringsAsFactors = FALSE),
    stratum_spec = list(variable = "condition", levels = "Mock"),
    dataset_id   = "coex_test"
  )

  res <- estimate_singlecellggm(bundle, n_iter = 100L, subsample = 20L,
                                 coex_cutoff = 10L, pcor_cutoff = 0.0,
                                 keep_negative = TRUE)
  et <- res[["Mock"]]$edge_table

  # All retained edges must satisfy coex_cells >= 10
  if (nrow(et) > 0L) {
    expect_true(all(et$coex_cells >= 10L),
                label = "all edges meet coex_cutoff = 10")
  }

  # Pair (gene1, gene2) specifically excluded (5 co-detected cells < 10)
  g1 <- gene_ids[1L]; g2 <- gene_ids[2L]
  pair_present <- nrow(et) > 0L &&
    any((et$gene_id_A == g1 & et$gene_id_B == g2) |
        (et$gene_id_A == g2 & et$gene_id_B == g1))
  expect_false(pair_present, label = "low-coex pair excluded from edge_table")
})

test_that("singlecellggm: keep_negative = FALSE excludes negative-weight edges", {
  bundle <- .make_ggm_bundle()
  res    <- estimate_singlecellggm(bundle, n_iter = 10L, subsample = 15L,
                                    pcor_cutoff = 0.0, keep_negative = FALSE)
  for (nr in res) {
    et <- nr$edge_table
    if (nrow(et) > 0L) {
      expect_true(all(et$weight >= 0),
                  label = paste("all weights >= 0 with keep_negative=FALSE for", nr$stratum_id))
    }
  }
})

test_that("singlecellggm: keep_negative = TRUE retains negative-weight edges", {
  bundle <- .make_ggm_bundle()
  res    <- estimate_singlecellggm(bundle, n_iter = 10L, subsample = 15L,
                                    pcor_cutoff = 0.0, keep_negative = TRUE)
  all_weights <- unlist(lapply(res, function(nr) nr$edge_table$weight))
  expect_true(any(all_weights < 0),
              label = "negative weights present with keep_negative=TRUE and pcor_cutoff=0")
})

test_that("singlecellggm: all gene IDs in edge_table are in bundle$gene_meta$gene_id", {
  bundle    <- .make_ggm_bundle()
  res       <- estimate_singlecellggm(bundle, n_iter = 3L, subsample = 10L,
                                       pcor_cutoff = 0.0, keep_negative = TRUE)
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

test_that("singlecellggm: p > subsample path runs correctly (real iteration loop)", {
  # n_genes = 40, subsample = 10 → 40 > 10 → real sampling loop (not guard path)
  bundle <- .make_ggm_bundle(n_genes = 40L)
  res    <- estimate_singlecellggm(bundle, n_iter = 5L, subsample = 10L,
                                    pcor_cutoff = 0.0, keep_negative = TRUE)
  for (nr in res) {
    expect_true(.is_network_result(nr),
                label = paste("valid NetworkResult for", nr$stratum_id))
    expect_equal(nr$params$subsample, 10L)
    expect_true(nr$params$n_genes > 10L,
                label = "p > subsample confirmed for real loop path")
  }
})

test_that("singlecellggm: does not load corpcor as a side-effect", {
  bundle    <- .make_ggm_bundle(n_genes = 40L)
  ns_before <- loadedNamespaces()
  estimate_singlecellggm(bundle, n_iter = 2L, subsample = 10L)
  ns_after  <- loadedNamespaces()
  new_ns    <- setdiff(ns_after, ns_before)
  expect_false("corpcor" %in% new_ns,
               label = "corpcor not loaded by estimate_singlecellggm")
})
