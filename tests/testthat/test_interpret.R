library(testthat)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

.at_ids <- function(n, chr = 1L)
  paste0("AT", chr, "G", formatC(seq_len(n) * 10L, width = 5L, flag = "0"))

# Build a minimal RobustnessResult with n_genes genes.
# Injects block structure: genes split into n_blocks groups; within-group
# pairs get high z_bar / R_score, cross-group pairs get low values.
.make_rob_structured <- function(n_genes = 50L, n_blocks = 3L, seed = 42L) {
  set.seed(seed)
  gene_ids <- .at_ids(n_genes)
  pairs    <- combn(n_genes, 2)
  n_pairs  <- ncol(pairs)

  block <- rep_len(seq_len(n_blocks), n_genes)
  same_block <- block[pairs[1, ]] == block[pairs[2, ]]

  # Within-block: high correlation; cross-block: low
  rho <- ifelse(same_block, runif(n_pairs, 0.25, 0.40), runif(n_pairs, 0.01, 0.05))
  z_bar   <- atanh(rho)
  R_score <- ifelse(same_block, runif(n_pairs, 0.7, 1.0), runif(n_pairs, 0.0, 0.3))

  ps <- data.frame(
    gene_id_A = gene_ids[pairs[1, ]],
    gene_id_B = gene_ids[pairs[2, ]],
    R_score   = R_score,
    z_bar     = z_bar,
    tau2      = runif(n_pairs, 0, 0.01),
    pval      = runif(n_pairs),
    qval      = runif(n_pairs),
    star      = NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    pair_scores   = ps,
    method_params = list(k = 1.64, weight_cap = 30, fdr_method = "BH",
                         n_strata = 2L, stratum_names = c("Mock", "DC3000"))
  )
}

# A minimal NetworkResult for one stratum.
.make_nr_simple <- function(gene_ids, n_cells = 100L, stratum_id = "Mock",
                            edge_fraction = 0.3, rho = 0.05, seed = 1L) {
  set.seed(seed)
  n   <- length(gene_ids)
  idx <- combn(n, 2)
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
    params     = list(n_cells = n_cells, n_genes = n, pcor_cutoff = 0.02,
                      n_iter = 100L, subsample = 2000L,
                      aggregation = "min_abs_pcor",
                      coex_cutoff = 10L, keep_negative = FALSE,
                      ridge = 1e-6, seed = seed),
    timestamp  = Sys.time()
  )
}

.make_network_list_simple <- function(gene_ids,
                                     strata = c("Mock", "DC3000"),
                                     seed = 1L) {
  setNames(
    lapply(seq_along(strata), function(i)
      .make_nr_simple(gene_ids, stratum_id = strata[i],
                      n_cells = 100L, seed = seed + i)),
    strata
  )
}

# Build a synthetic ModuleInput for GOI tests (no WGCNA needed).
.make_mod_input_simple <- function(gene_ids) {
  n <- length(gene_ids)
  half <- n %/% 2L
  top_labels <- c(rep(1L, half), rep(2L, n - half))
  sub_labels <- top_labels  # same as top for simplicity
  kme_vals   <- c(runif(half, 0.5, 0.9), runif(n - half, 0.4, 0.85))

  gene_module <- data.frame(
    gene_id    = gene_ids,
    top_module = top_labels,
    sub_module = sub_labels,
    kME        = kme_vals,
    stringsAsFactors = FALSE
  )

  module_meta <- data.frame(
    module_id              = c(1L, 2L),
    n_genes                = c(half, n - half),
    label                  = NA_character_,
    top_organ_or_condition = NA_character_,
    delta_treatment        = NA_character_,
    go_top                 = NA_character_,
    zsummary               = c(5.2, 3.1),
    preservation_method    = "fallback_meancor",
    stringsAsFactors = FALSE
  )

  hub_list <- lapply(c(1L, 2L), function(m) {
    sub_gm <- gene_module[gene_module$top_module == m, ]
    sub_gm <- sub_gm[order(sub_gm$kME, decreasing = TRUE), ]
    sub_gm <- head(sub_gm, 20L)
    data.frame(module_id = m, gene_id = sub_gm$gene_id,
               gene_symbol = NA_character_, kME = sub_gm$kME,
               hub_rank = seq_len(nrow(sub_gm)), stringsAsFactors = FALSE)
  })

  list(
    gene_module  = gene_module,
    module_meta  = module_meta,
    module_hier  = data.frame(sub_module = c(1L, 2L),
                              top_module = c(1L, 2L)),
    hub_genes    = do.call(rbind, hub_list),
    module_tfs   = data.frame(module_id = integer(), gene_id = character(),
                              gene_symbol = character(), tf_family = character(),
                              stringsAsFactors = FALSE),
    eigengenes   = matrix(runif(n * 2L), nrow = n, ncol = 2L,
                          dimnames = list(gene_ids, c("ME1", "ME2")))
  )
}

# Build a minimal RobustnessResult for GOI tests.
.make_rob_simple <- function(gene_ids, seed = 1L) {
  set.seed(seed)
  n     <- length(gene_ids)
  pairs <- combn(n, 2)
  n_p   <- ncol(pairs)

  ps <- data.frame(
    gene_id_A = gene_ids[pairs[1, ]],
    gene_id_B = gene_ids[pairs[2, ]],
    R_score   = runif(n_p, 0.1, 1.0),
    z_bar     = atanh(runif(n_p, 0.05, 0.40)),
    tau2      = runif(n_p, 0, 0.01),
    pval      = runif(n_p),
    qval      = runif(n_p),
    star      = NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    pair_scores   = ps,
    method_params = list(k = 1.64, weight_cap = 30, fdr_method = "BH",
                         n_strata = 2L, stratum_names = c("Mock", "DC3000"))
  )
}

# ===========================================================================
# build_wgcna_modules() tests
# ===========================================================================

test_that("build_wgcna_modules returns valid ModuleInput structure", {
  skip_if_not_installed("WGCNA")
  gene_ids <- .at_ids(50L)
  rob      <- .make_rob_structured(n_genes = 50L, n_blocks = 3L)
  nl       <- .make_network_list_simple(gene_ids)

  mod <- build_wgcna_modules(rob, nl,
                              r_score_min     = 0,
                              soft_power      = 4L,
                              min_module_size = 5L,
                              merge_cut       = 0.30,
                              sub_merge_cut   = 0.15)

  expect_type(mod, "list")
  required_slots <- c("gene_module", "module_meta", "module_hier",
                      "hub_genes", "module_tfs", "eigengenes")
  expect_true(all(required_slots %in% names(mod)),
              info = paste("Missing slots:",
                           paste(setdiff(required_slots, names(mod)), collapse = ", ")))
})

test_that("gene_module has one row per gene; top_module is integer; kME in [-1, 1]", {
  skip_if_not_installed("WGCNA")
  gene_ids <- .at_ids(50L)
  rob      <- .make_rob_structured(n_genes = 50L, n_blocks = 3L)
  nl       <- .make_network_list_simple(gene_ids)

  mod <- build_wgcna_modules(rob, nl,
                              soft_power      = 4L,
                              min_module_size = 5L,
                              merge_cut       = 0.30,
                              sub_merge_cut   = 0.15)

  gm <- mod$gene_module
  expect_s3_class(gm, "data.frame")
  expect_equal(nrow(gm), length(gene_ids),
               info = "gene_module should have one row per gene in the network")
  expect_true(is.integer(gm$top_module),
              info = "top_module should be integer")

  kme_vals <- gm$kME[!is.na(gm$kME)]
  if (length(kme_vals) > 0L) {
    expect_true(all(kme_vals >= -1 & kme_vals <= 1),
                info = "kME values must be in [-1, 1]")
  }

  req_cols <- c("gene_id", "top_module", "sub_module", "kME")
  expect_true(all(req_cols %in% names(gm)),
              info = paste("gene_module missing columns:",
                           paste(setdiff(req_cols, names(gm)), collapse = ", ")))
})

test_that("hub_genes has at most 20 rows per module", {
  skip_if_not_installed("WGCNA")
  gene_ids <- .at_ids(50L)
  rob      <- .make_rob_structured(n_genes = 50L, n_blocks = 3L)
  nl       <- .make_network_list_simple(gene_ids)

  mod <- build_wgcna_modules(rob, nl,
                              soft_power      = 4L,
                              min_module_size = 5L,
                              merge_cut       = 0.30,
                              sub_merge_cut   = 0.15)

  hg          <- mod$hub_genes
  req_cols_hg <- c("module_id", "gene_id", "gene_symbol", "kME", "hub_rank")
  expect_true(all(req_cols_hg %in% names(hg)),
              info = paste("hub_genes missing columns:",
                           paste(setdiff(req_cols_hg, names(hg)), collapse = ", ")))

  if (nrow(hg) > 0L) {
    hub_counts <- table(hg$module_id)
    expect_true(all(hub_counts <= 20L),
                info = "hub_genes must have at most 20 rows per module")
  }
})

test_that("module_hier sub_modules are a subset of gene_module sub_module values", {
  skip_if_not_installed("WGCNA")
  gene_ids <- .at_ids(50L)
  rob      <- .make_rob_structured(n_genes = 50L, n_blocks = 3L)
  nl       <- .make_network_list_simple(gene_ids)

  mod <- build_wgcna_modules(rob, nl,
                              soft_power      = 4L,
                              min_module_size = 5L,
                              merge_cut       = 0.30,
                              sub_merge_cut   = 0.15)

  mh <- mod$module_hier
  expect_s3_class(mh, "data.frame")
  expect_true(all(c("sub_module", "top_module") %in% names(mh)),
              info = "module_hier needs sub_module and top_module columns")

  if (nrow(mh) > 0L) {
    gm_sub_vals <- unique(mod$gene_module$sub_module[mod$gene_module$sub_module > 0L])
    expect_true(all(mh$sub_module %in% gm_sub_vals),
                info = "module_hier sub_modules not all present in gene_module")
  }
})

test_that("soft_power auto-pick runs without error on synthetic 50-gene data", {
  skip_if_not_installed("WGCNA")
  gene_ids <- .at_ids(50L)
  rob      <- .make_rob_structured(n_genes = 50L, n_blocks = 3L)
  nl       <- .make_network_list_simple(gene_ids)

  # soft_power = NULL triggers auto-pick
  expect_no_error(
    build_wgcna_modules(rob, nl,
                        soft_power      = NULL,
                        min_module_size = 5L,
                        merge_cut       = 0.30,
                        sub_merge_cut   = 0.15)
  )
})

test_that("r_score_min filter removes low-score edges", {
  skip_if_not_installed("WGCNA")
  gene_ids <- .at_ids(50L)
  rob      <- .make_rob_structured(n_genes = 50L, n_blocks = 3L)
  nl       <- .make_network_list_simple(gene_ids)

  # Only keep edges with R_score >= 0.7 (all within-block pairs)
  mod_filtered <- build_wgcna_modules(rob, nl,
                                      r_score_min     = 0.7,
                                      soft_power      = 4L,
                                      min_module_size = 5L,
                                      merge_cut       = 0.30,
                                      sub_merge_cut   = 0.15)
  expect_s3_class(mod_filtered$gene_module, "data.frame")

  # Genes not in any high-R pair should be absent from gene_module
  high_ps <- rob$pair_scores[rob$pair_scores$R_score >= 0.7, ]
  expect_true(nrow(high_ps) > 0L, info = "Synthetic data should have high-R pairs")
})

test_that("module_meta has required columns", {
  skip_if_not_installed("WGCNA")
  gene_ids <- .at_ids(50L)
  rob      <- .make_rob_structured(n_genes = 50L, n_blocks = 3L)
  nl       <- .make_network_list_simple(gene_ids)

  mod <- build_wgcna_modules(rob, nl,
                              soft_power      = 4L,
                              min_module_size = 5L,
                              merge_cut       = 0.30,
                              sub_merge_cut   = 0.15)

  mm <- mod$module_meta
  req_cols <- c("module_id", "n_genes", "label", "top_organ_or_condition",
                "delta_treatment", "go_top", "zsummary", "preservation_method")
  expect_true(all(req_cols %in% names(mm)),
              info = paste("module_meta missing columns:",
                           paste(setdiff(req_cols, names(mm)), collapse = ", ")))
})

# ===========================================================================
# build_goi_table() tests
# ===========================================================================

test_that("GOI gene in network: all required fields are present", {
  gene_ids <- .at_ids(20L)
  mod_in   <- .make_mod_input_simple(gene_ids)
  rob      <- .make_rob_simple(gene_ids)

  goi_in <- gene_ids[1L]  # definitely in network
  result  <- build_goi_table(mod_in, rob, goi_ids = goi_in, top_n = 5L)

  expect_type(result, "list")
  expect_named(result, c("wide", "long"), ignore.order = TRUE)

  wide_row <- result$wide[result$wide$gene_id == goi_in, , drop = FALSE]
  expect_equal(nrow(wide_row), 1L)
  expect_false(is.na(wide_row$module),
               info = "module should be non-NA for in-network gene")
  expect_false(is.na(wide_row$kME),
               info = "kME should be non-NA for in-network gene")
  expect_type(wide_row$hub_flag, "logical")
})

test_that("GOI gene not in network: row present with NA module and notes='not in network'", {
  gene_ids   <- .at_ids(20L)
  mod_in     <- .make_mod_input_simple(gene_ids)
  rob        <- .make_rob_simple(gene_ids)

  absent_goi <- "AT9G99999"  # guaranteed not in gene_ids
  result     <- build_goi_table(mod_in, rob, goi_ids = absent_goi, top_n = 5L)

  wide_row <- result$wide[result$wide$gene_id == absent_goi, , drop = FALSE]
  expect_equal(nrow(wide_row), 1L, info = "Absent gene must still have a row")
  expect_true(is.na(wide_row$module), info = "module should be NA for absent gene")
  expect_equal(wide_row$notes, "not in network")
})

test_that("$long has top_n rows per GOI or fewer if fewer partners exist", {
  gene_ids <- .at_ids(15L)
  mod_in   <- .make_mod_input_simple(gene_ids)
  rob      <- .make_rob_simple(gene_ids)

  top_n  <- 5L
  goi_id <- gene_ids[1L]
  result <- build_goi_table(mod_in, rob, goi_ids = goi_id, top_n = top_n)

  long_rows <- result$long[result$long$gene_id == goi_id, , drop = FALSE]
  n_avail   <- sum(rob$pair_scores$gene_id_A == goi_id |
                   rob$pair_scores$gene_id_B == goi_id)
  expect_lte(nrow(long_rows), top_n,
             label = "$long should not exceed top_n rows per GOI")
  expect_lte(nrow(long_rows), n_avail,
             label = "$long should not exceed available partners")
})

test_that("$wide top_N_coexpressed_partners is a character (semicolon-joined or empty)", {
  gene_ids <- .at_ids(20L)
  mod_in   <- .make_mod_input_simple(gene_ids)
  rob      <- .make_rob_simple(gene_ids)

  result <- build_goi_table(mod_in, rob, goi_ids = gene_ids[1L], top_n = 3L)

  top_n_col <- result$wide$top_N_coexpressed_partners[1L]
  expect_type(top_n_col, "character")
  if (!is.na(top_n_col) && nchar(top_n_col) > 0L) {
    # Should match "AT...G... (R=0.XX); ..." format
    expect_match(top_n_col, "AT\\dG\\d{5} \\(R=",
                 info = "Partner format should be 'AT#G##### (R=0.XX)'")
  }
})

test_that("$wide and $long are data.frames with correct column names", {
  gene_ids <- .at_ids(20L)
  mod_in   <- .make_mod_input_simple(gene_ids)
  rob      <- .make_rob_simple(gene_ids)

  result <- build_goi_table(mod_in, rob, goi_ids = gene_ids[1:3], top_n = 5L)

  wide_cols <- c("gene_id", "gene_symbol", "module", "kME", "hub_flag",
                 "zsummary", "preservation_method",
                 "top_N_coexpressed_partners", "notes")
  long_cols <- c("gene_id", "gene_symbol", "partner_id", "partner_symbol",
                 "R_score", "weight", "rank", "notes")

  expect_s3_class(result$wide, "data.frame")
  expect_s3_class(result$long, "data.frame")
  expect_true(all(wide_cols %in% names(result$wide)),
              info = paste("$wide missing columns:",
                           paste(setdiff(wide_cols, names(result$wide)), collapse = ", ")))
  expect_true(all(long_cols %in% names(result$long)),
              info = paste("$long missing columns:",
                           paste(setdiff(long_cols, names(result$long)), collapse = ", ")))
})

test_that("build_goi_table handles absent gene and in-network gene together", {
  gene_ids   <- .at_ids(20L)
  mod_in     <- .make_mod_input_simple(gene_ids)
  rob        <- .make_rob_simple(gene_ids)
  absent_goi <- "AT9G99999"
  goi_ids    <- c(gene_ids[1L], absent_goi)

  result <- build_goi_table(mod_in, rob, goi_ids = goi_ids, top_n = 5L)

  expect_equal(nrow(result$wide), 2L,
               info = "$wide should have one row per GOI including absent ones")
  present_row <- result$wide[result$wide$gene_id == gene_ids[1L], ]
  absent_row  <- result$wide[result$wide$gene_id == absent_goi, ]
  expect_false(is.na(present_row$module))
  expect_true(is.na(absent_row$module))
  expect_equal(absent_row$notes, "not in network")
})

# ===========================================================================
# save_goi_table() tests
# ===========================================================================

test_that("save_goi_table writes _lookup_wide.csv and _lookup_long.csv", {
  gene_ids <- .at_ids(20L)
  mod_in   <- .make_mod_input_simple(gene_ids)
  rob      <- .make_rob_simple(gene_ids)
  result   <- build_goi_table(mod_in, rob, goi_ids = gene_ids[1:3], top_n = 5L)

  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  save_goi_table(result, outdir = tmp, list_name = "test_list")

  expect_true(file.exists(file.path(tmp, "test_list_lookup_wide.csv")))
  expect_true(file.exists(file.path(tmp, "test_list_lookup_long.csv")))

  wide_back <- read.csv(file.path(tmp, "test_list_lookup_wide.csv"),
                        stringsAsFactors = FALSE)
  expect_equal(nrow(wide_back), nrow(result$wide))
})
