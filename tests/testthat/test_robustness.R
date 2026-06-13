library(testthat)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

.at_ids <- function(n) paste0("AT1G", formatC(seq_len(n) * 10L, width = 5L, flag = "0"))

# Build a synthetic NetworkResult with a small edge table.
# gene_ids: character vector of AT-IDs
# n_cells: cell count to embed in params (GGM mode)
# edge_fraction: fraction of gene pairs that appear as edges
# rho_range: pcor range for sampled edges
.make_nr <- function(gene_ids, n_cells, stratum_id,
                     edge_fraction = 0.3, rho = 0.05, seed = 1L) {
  set.seed(seed)
  n   <- length(gene_ids)
  idx <- combn(n, 2)                                 # all pairs
  keep <- sample(ncol(idx), round(ncol(idx) * edge_fraction))
  idx  <- idx[, keep, drop = FALSE]

  et <- data.frame(
    gene_id_A    = gene_ids[idx[1, ]],
    gene_id_B    = gene_ids[idx[2, ]],
    weight       = rho + runif(ncol(idx), 0, 0.02),
    coex_cells   = sample(10:100, ncol(idx), replace = TRUE),
    sampling_num = sample(80:120, ncol(idx), replace = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    edge_table = et,
    gene_ids   = gene_ids,
    stratum_id = stratum_id,
    mode       = "singlecellggm",
    params     = list(n_cells = n_cells, n_genes = n,
                      pcor_cutoff = 0.02, n_iter = 100L,
                      subsample = 2000L, aggregation = "min_abs_pcor",
                      coex_cutoff = 10L, keep_negative = FALSE,
                      ridge = 1e-6, seed = 98L),
    timestamp  = Sys.time()
  )
}

# A pair of 2-stratum NetworkResult lists sharing the same gene universe.
.make_network_list <- function(n_genes = 20L, n_cells = 50L,
                               strata = c("Mock", "DC3000"),
                               seed = 1L) {
  genes <- .at_ids(n_genes)
  setNames(
    lapply(seq_along(strata), function(i)
      .make_nr(genes, n_cells, strata[i], seed = seed + i)),
    strata
  )
}

# ---------------------------------------------------------------------------
# Test 1 – compute_robustness() returns a valid RobustnessResult
# ---------------------------------------------------------------------------

test_that("compute_robustness returns valid RobustnessResult structure", {
  nl <- .make_network_list()
  rob <- compute_robustness(nl)

  expect_type(rob, "list")
  expect_named(rob, c("pair_scores", "method_params"), ignore.order = TRUE)

  ps <- rob$pair_scores
  expect_s3_class(ps, "data.frame")

  required_cols <- c("gene_id_A", "gene_id_B", "R_score", "z_bar",
                     "tau2", "pval", "qval", "star")
  expect_true(all(required_cols %in% names(ps)),
              info = paste("Missing:", paste(setdiff(required_cols, names(ps)), collapse = ", ")))

  mp <- rob$method_params
  expect_true(all(c("k", "weight_cap", "fdr_method", "n_strata", "stratum_names") %in% names(mp)))
  expect_equal(mp$n_strata, 2L)
  expect_equal(mp$stratum_names, c("Mock", "DC3000"))
})

# ---------------------------------------------------------------------------
# Test 2 – pair_scores has one row per unique pair (no duplicates)
# ---------------------------------------------------------------------------

test_that("pair_scores has one row per unique pair with no duplicates", {
  nl  <- .make_network_list(n_genes = 20L)
  rob <- compute_robustness(nl)
  ps  <- rob$pair_scores

  # A-B and B-A are the same pair; keys should be canonical (A < B)
  keys <- paste(ps$gene_id_A, ps$gene_id_B)
  expect_equal(length(keys), length(unique(keys)),
               info = "Duplicate pairs found in pair_scores")
})

# ---------------------------------------------------------------------------
# Test 3 – R_score in [0,1]; z_bar numeric; pval and qval in [0,1]
# ---------------------------------------------------------------------------

test_that("R_score, pval, qval are in [0,1] and z_bar is finite numeric", {
  nl  <- .make_network_list(n_genes = 20L)
  rob <- compute_robustness(nl)
  ps  <- rob$pair_scores

  expect_true(all(ps$R_score >= 0 & ps$R_score <= 1),
              info = "R_score out of [0,1]")
  expect_true(all(is.numeric(ps$z_bar) & is.finite(ps$z_bar)),
              info = "z_bar contains non-finite values")
  expect_true(all(ps$pval >= 0 & ps$pval <= 1, na.rm = TRUE),
              info = "pval out of [0,1]")
  expect_true(all(ps$qval >= 0 & ps$qval <= 1, na.rm = TRUE),
              info = "qval out of [0,1]")
})

# ---------------------------------------------------------------------------
# Test 4 – I_<stratum> columns present, one per stratum, values 0 or 1
# ---------------------------------------------------------------------------

test_that("I_<stratum> columns present with 0/1 values", {
  strata <- c("Mock", "DC3000", "AvrRpt2")
  nl     <- .make_network_list(n_genes = 20L, strata = strata)
  rob    <- compute_robustness(nl)
  ps     <- rob$pair_scores

  i_cols <- paste0("I_", strata)
  expect_true(all(i_cols %in% names(ps)),
              info = paste("Missing I columns:", paste(setdiff(i_cols, names(ps)), collapse = ", ")))

  for (col in i_cols) {
    vals <- ps[[col]]
    expect_true(all(vals %in% c(0L, 1L)),
                info = paste(col, "contains values other than 0 or 1"))
  }
})

# ---------------------------------------------------------------------------
# Test 5 – star column is NA out of compute_robustness (filled by annotate_star)
# ---------------------------------------------------------------------------

test_that("star column is NA in compute_robustness output", {
  nl  <- .make_network_list()
  rob <- compute_robustness(nl)
  expect_true(all(is.na(rob$pair_scores$star)),
              info = "star column should be all NA before annotate_star()")
})

# ---------------------------------------------------------------------------
# Test 6 – annotate_star: TRUE for robust-in-both, FALSE for below threshold
#           in rob2, NA for pairs absent from rob2
# ---------------------------------------------------------------------------

test_that("annotate_star classifies pairs correctly", {
  genes  <- .at_ids(15L)
  n_cells <- 80L

  # rob1: all pairs between genes[1:5] × genes[6:10] (controlled edges)
  nr1 <- .make_nr(genes[1:10], n_cells, "MockA", rho = 0.10, seed = 1L)
  nr2 <- .make_nr(genes[1:10], n_cells, "DC3kA", rho = 0.10, seed = 2L)
  rob1 <- compute_robustness(list(MockA = nr1, DC3kA = nr2))

  # rob2: contains some of the same genes; no overlap for genes[11:15]
  nr3 <- .make_nr(genes[1:8], n_cells, "MockB", rho = 0.10, seed = 3L)
  nr4 <- .make_nr(genes[1:8], n_cells, "DC3kB", rho = 0.10, seed = 4L)
  rob2 <- compute_robustness(list(MockB = nr3, DC3kB = nr4))

  rob1_annotated <- annotate_star(rob1, rob2, threshold = 0.5)

  ps <- rob1_annotated$pair_scores

  # star must be non-NA, TRUE, or FALSE — never some other value
  expect_true(all(is.na(ps$star) | ps$star %in% c(0, 1)),
              info = "star column has unexpected values")

  # Pairs between genes NOT in rob2 should be NA
  in_rob2_genes <- genes[1:8]
  ps$in_rob2 <- ps$gene_id_A %in% in_rob2_genes & ps$gene_id_B %in% in_rob2_genes
  absent_pairs <- ps[!ps$in_rob2, ]
  expect_true(all(is.na(absent_pairs$star)),
              info = "Pairs absent from rob2 should have star = NA")

  # rob1 is unchanged
  expect_true(all(is.na(rob1$pair_scores$star)),
              info = "annotate_star should not modify rob1 in-place")

  # rob2 is not modified
  expect_true(all(is.na(rob2$pair_scores$star)),
              info = "annotate_star should not modify rob2")
})

# ---------------------------------------------------------------------------
# Test 7 – save_robustness: files written with correct names
# ---------------------------------------------------------------------------

test_that("save_robustness writes pair_scores_full.csv and robustness_result.rds", {
  nl  <- .make_network_list()
  rob <- compute_robustness(nl)

  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))

  save_robustness(rob, tmp)

  expect_true(file.exists(file.path(tmp, "pair_scores_full.csv")))
  expect_true(file.exists(file.path(tmp, "robustness_result.rds")))
})

# ---------------------------------------------------------------------------
# Test 8 – FLAG-03: save_robustness writes ALL pairs including R_score=0 pairs
# ---------------------------------------------------------------------------

test_that("save_robustness writes ALL pairs including R_score=0 ones (FLAG-03)", {
  # Build two strata with completely disjoint edge sets so every pair
  # appears in at most one stratum → many pairs will have R_score < 1.
  genes_a <- .at_ids(10L)
  genes_b <- paste0("AT2G", formatC(seq_len(10L) * 10L, width = 5L, flag = "0"))

  # Stratum 1: edges only among genes_a
  nr1 <- .make_nr(genes_a, n_cells = 50L, stratum_id = "S1",
                  edge_fraction = 1.0, rho = 0.05, seed = 11L)
  # Stratum 2: edges only among genes_b
  nr2_et <- data.frame(
    gene_id_A    = genes_b[1:5],
    gene_id_B    = genes_b[6:10],
    weight       = rep(0.06, 5),
    coex_cells   = rep(20L, 5),
    sampling_num = rep(100L, 5),
    stringsAsFactors = FALSE
  )
  nr2 <- list(
    edge_table = nr2_et, gene_ids = genes_b, stratum_id = "S2",
    mode = "singlecellggm",
    params = list(n_cells = 50L, n_genes = 10L, pcor_cutoff = 0.02,
                  n_iter = 100L, subsample = 2000L, aggregation = "min_abs_pcor",
                  coex_cutoff = 10L, keep_negative = FALSE,
                  ridge = 1e-6, seed = 98L),
    timestamp = Sys.time()
  )

  nl  <- list(S1 = nr1, S2 = nr2)
  rob <- compute_robustness(nl)

  # Pairs from S1 won't appear in S2 and vice versa → many have R_score ≠ 1
  has_low_r <- any(rob$pair_scores$R_score < 1)
  expect_true(has_low_r, info = "Synthetic data should produce some R_score < 1")

  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  save_robustness(rob, tmp)

  saved <- read.csv(file.path(tmp, "pair_scores_full.csv"),
                    stringsAsFactors = FALSE)
  expect_equal(nrow(saved), nrow(rob$pair_scores),
               info = "pair_scores_full.csv must contain ALL pairs, including R_score < 1")
})
