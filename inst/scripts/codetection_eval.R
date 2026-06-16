#!/usr/bin/env Rscript
# Reusable single-cell co-detection evaluation engine
# Usage: source("inst/scripts/codetection_eval.R")
# Main export: score_gene_group()
# No data.table; base-R + Matrix sparse ops throughout.

suppressPackageStartupMessages(library(Matrix))

# ── Helper: compute detection rate for all genes ─────────────────────────────
compute_det_rates <- function(counts_mat) {
  # counts_mat: sparse genes × cells (dgCMatrix or similar)
  # Returns named numeric vector (gene → detection_rate)
  n_cells <- ncol(counts_mat)
  bool_mat <- counts_mat > 0
  dr <- Matrix::rowSums(bool_mat) / n_cells
  setNames(as.numeric(dr), rownames(counts_mat))
}

# ── Helper: safe quantile breaks (handles ties/duplicates) ──────────────────
safe_breaks <- function(x, n_bins = 10L) {
  qs <- quantile(x, probs = seq(0, 1, length.out = n_bins + 1L), names = FALSE)
  qs[1]              <- qs[1] - 1e-9
  qs[length(qs)]     <- qs[length(qs)] + 1e-9
  # Remove duplicates while preserving order
  qs <- unique(qs)
  qs
}

# ── Helper: build expression-frequency decile pools ─────────────────────────
build_decile_pools <- function(det_rate_all, exclude_genes = character(0)) {
  # Returns list indexed by bin number; each element = gene IDs in that bin
  brks        <- safe_breaks(det_rate_all)
  all_bins    <- cut(det_rate_all, breaks = brks, labels = FALSE, include.lowest = TRUE)
  names(all_bins) <- names(det_rate_all)
  n_bins      <- length(brks) - 1L
  excl_set    <- setNames(rep(TRUE, length(exclude_genes)), exclude_genes)
  lapply(seq_len(n_bins), function(d) {
    cands    <- names(all_bins)[!is.na(all_bins) & all_bins == d]
    in_excl  <- cands %in% names(excl_set)
    if (sum(!in_excl) > 0L) cands[!in_excl] else cands
  })
}

# ── Core function ─────────────────────────────────────────────────────────────
#
# score_gene_group()
#   counts_mat   : sparse Matrix, genes × cells (raw counts), rownames = RNA gene IDs
#   lognorm_mat  : sparse Matrix, genes × cells (log-normalized), same rownames
#   det_rate_all : named numeric from compute_det_rates(counts_mat)
#   gene_group   : character vector of rownames(counts_mat) to evaluate (pre-resolved)
#   group_id     : label for logging
#   n_pair_cap   : max pairs sampled (default 5000)
#   seed         : RNG seed
#   timeout_secs : per-group wall-clock limit (seconds); returns partial if exceeded
#   log_fn       : message function
#
# Returns: list with summary stats + pair-level vectors for further analysis.
# Returns NULL if group has < 2 valid genes.

score_gene_group <- function(
    counts_mat,
    lognorm_mat,
    det_rate_all,
    gene_group,
    group_id     = "group",
    n_pair_cap   = 5000L,
    seed         = 98L,
    timeout_secs = 600L,
    log_fn       = message
) {
  set.seed(seed)
  t_start <- proc.time()[["elapsed"]]

  # ── Restrict to genes present in matrix ─────────────────────────────────
  avail <- gene_group[gene_group %in% rownames(counts_mat)]
  n_miss <- length(gene_group) - length(avail)
  if (n_miss > 0)
    log_fn(sprintf("  [%s] %d/%d genes missing from matrix", group_id, n_miss, length(gene_group)))
  if (length(avail) < 2L) {
    log_fn(sprintf("  [%s] < 2 genes available — skipping", group_id))
    return(NULL)
  }

  n_genes  <- length(avail)
  n_cells  <- ncol(counts_mat)

  # ── Sample pairs ─────────────────────────────────────────────────────────
  # Use choose() to check total without enumerating — combn() on large n is fatal.
  n_total <- choose(n_genes, 2L)
  capped  <- n_total > n_pair_cap
  if (capped) {
    # Direct random sampling of pairs (no enumeration)
    i1 <- sample.int(n_genes, n_pair_cap, replace = TRUE)
    i2 <- sample.int(n_genes - 1L, n_pair_cap, replace = TRUE)
    i2[i2 >= i1] <- i2[i2 >= i1] + 1L   # ensure i2 ≠ i1
    all_i <- rbind(pmin(i1, i2), pmax(i1, i2))
    log_fn(sprintf("  [%s] %d genes → capped to %d/%d pairs (direct sampling)",
                   group_id, n_genes, n_pair_cap, n_total))
  } else {
    all_i <- combn(n_genes, 2L)           # safe when n_genes is small
  }
  n_pairs <- ncol(all_i)
  g1 <- avail[all_i[1L, ]]
  g2 <- avail[all_i[2L, ]]

  # ── Pre-extract unique group genes for efficiency ─────────────────────────
  unique_grp <- unique(c(g1, g2))   # at most n_genes (bounded)
  n_unique   <- length(unique_grp)

  # For groups ≤ 300 unique genes: dense extraction + rank-transform once
  use_dense <- n_unique <= 300L
  if (use_dense) {
    grp_ln   <- as.matrix(lognorm_mat[unique_grp, , drop = FALSE])  # n_unique × n_cells
    grp_bool <- as.matrix(counts_mat[unique_grp, , drop = FALSE] > 0)
    gu_idx   <- setNames(seq_len(n_unique), unique_grp)
    # rank-transform for fast Spearman (rank once, compute Pearson of ranks)
    grp_rank <- matrix(0, nrow = n_unique, ncol = n_cells)
    for (gi in seq_len(n_unique))
      grp_rank[gi, ] <- rank(grp_ln[gi, ], ties.method = "average")
    # Precompute rank means & sds for Pearson
    r_means <- rowMeans(grp_rank)
    r_sds   <- apply(grp_rank, 1, sd)
    r_sds[r_sds == 0] <- NA_real_
  }

  # ── Within-group co-detection + Spearman ─────────────────────────────────
  codet_within <- numeric(n_pairs)
  spear_within <- numeric(n_pairs)
  timed_out <- FALSE

  for (k in seq_len(n_pairs)) {
    if (k %% 500L == 0L) {
      elapsed <- proc.time()[["elapsed"]] - t_start
      if (elapsed > timeout_secs) {
        log_fn(sprintf("  [%s] TIMED OUT at pair %d/%d (%.0fs)", group_id, k, n_pairs, elapsed))
        n_pairs    <- k - 1L
        g1         <- g1[seq_len(n_pairs)]
        g2         <- g2[seq_len(n_pairs)]
        codet_within <- codet_within[seq_len(n_pairs)]
        spear_within <- spear_within[seq_len(n_pairs)]
        timed_out  <- TRUE
        break
      }
    }

    if (use_dense) {
      i1 <- gu_idx[g1[k]]; i2 <- gu_idx[g2[k]]
      b1 <- grp_bool[i1, ]; b2 <- grp_bool[i2, ]
      codet_within[k] <- mean(b1 & b2)
      rc1 <- grp_rank[i1, ] - r_means[i1]
      rc2 <- grp_rank[i2, ] - r_means[i2]
      denom <- (n_cells - 1L) * r_sds[i1] * r_sds[i2]
      spear_within[k] <- if (!is.na(denom) && denom > 0) sum(rc1 * rc2) / denom else NA_real_
    } else {
      b1 <- as.logical(counts_mat[g1[k], ] > 0)
      b2 <- as.logical(counts_mat[g2[k], ] > 0)
      codet_within[k] <- mean(b1 & b2)
      x <- as.numeric(lognorm_mat[g1[k], ])
      y <- as.numeric(lognorm_mat[g2[k], ])
      spear_within[k] <- cor(x, y, method = "spearman")
    }
  }
  if (timed_out) n_pairs <- length(codet_within)

  # ── Build expression-frequency-matched null ───────────────────────────────
  brks_q     <- safe_breaks(det_rate_all)
  n_bins_q   <- length(brks_q) - 1L
  bin_pools  <- build_decile_pools(det_rate_all, exclude_genes = avail)
  # Assign each group gene to its bin
  g1_bin <- cut(det_rate_all[g1], breaks = brks_q, labels = FALSE, include.lowest = TRUE)
  g2_bin <- cut(det_rate_all[g2], breaks = brks_q, labels = FALSE, include.lowest = TRUE)

  # Draw null gene IDs for each pair
  mid_bin <- max(1L, round(n_bins_q / 2L))   # fallback if NA
  null_g1 <- character(n_pairs)
  null_g2 <- character(n_pairs)
  for (k in seq_len(n_pairs)) {
    d1 <- g1_bin[k]; d2 <- g2_bin[k]
    if (is.na(d1)) d1 <- mid_bin
    if (is.na(d2)) d2 <- mid_bin
    pool1 <- bin_pools[[d1]]; pool2 <- bin_pools[[d2]]
    if (length(pool1) == 0L) pool1 <- names(det_rate_all)
    if (length(pool2) == 0L) pool2 <- names(det_rate_all)
    null_g1[k] <- pool1[sample.int(length(pool1), 1L)]
    null_g2[k] <- pool2[sample.int(length(pool2), 1L)]
  }

  # Evaluate null pairs (always pair-by-pair to avoid large dense extraction)
  null_codet <- numeric(n_pairs)
  null_spear <- numeric(n_pairs)
  for (k in seq_len(n_pairs)) {
    b1 <- as.logical(counts_mat[null_g1[k], ] > 0)
    b2 <- as.logical(counts_mat[null_g2[k], ] > 0)
    null_codet[k] <- mean(b1 & b2)
    x <- as.numeric(lognorm_mat[null_g1[k], ])
    y <- as.numeric(lognorm_mat[null_g2[k], ])
    null_spear[k] <- cor(x, y, method = "spearman")
  }

  # ── Summarize ─────────────────────────────────────────────────────────────
  mc_w  <- mean(codet_within, na.rm = TRUE)
  mc_n  <- mean(null_codet,   na.rm = TRUE)
  gap_c <- mc_w - mc_n
  sd_nc <- sd(null_codet,     na.rm = TRUE)
  es_c  <- if (!is.na(sd_nc) && sd_nc > 0) gap_c / sd_nc else NA_real_

  ms_w  <- mean(spear_within, na.rm = TRUE)
  ms_n  <- mean(null_spear,   na.rm = TRUE)
  gap_s <- ms_w - ms_n
  sd_ns <- sd(null_spear,     na.rm = TRUE)
  es_s  <- if (!is.na(sd_ns) && sd_ns > 0) gap_s / sd_ns else NA_real_

  elapsed_total <- proc.time()[["elapsed"]] - t_start
  log_fn(sprintf("  [%s] Done: %d genes, %d pairs, codet_gap=%.4f (ES=%.2f), spear_gap=%.4f (ES=%.2f), %.0fs",
    group_id, n_genes, n_pairs, gap_c, ifelse(is.na(es_c), 0, es_c),
    gap_s, ifelse(is.na(es_s), 0, es_s), elapsed_total))

  list(
    group_id          = group_id,
    n_genes           = n_genes,
    n_pairs           = n_pairs,
    n_pairs_total     = n_total,
    capped            = capped,
    timed_out         = timed_out,
    mean_codet_within = mc_w,
    mean_codet_null   = mc_n,
    gap_codet         = gap_c,
    effect_size_codet = es_c,
    mean_spear_within = ms_w,
    mean_spear_null   = ms_n,
    gap_spear         = gap_s,
    effect_size_spear = es_s,
    codet_pairs       = codet_within,
    spear_pairs       = spear_within,
    null_codet_pairs  = null_codet,
    null_spear_pairs  = null_spear,
    g1_ids            = g1,
    g2_ids            = g2,
    null_g1_ids       = null_g1,
    null_g2_ids       = null_g2,
    elapsed_secs      = elapsed_total
  )
}

# ── Batch runner: score a named list of groups ─────────────────────────────
score_group_list <- function(
    group_list,     # named list: group_id -> character vector of RNA rownames
    counts_mat, lognorm_mat, det_rate_all,
    n_pair_cap = 5000L, seed = 98L, timeout_secs = 600L, log_fn = message
) {
  results <- vector("list", length(group_list))
  names(results) <- names(group_list)
  for (gid in names(group_list)) {
    log_fn(sprintf("Scoring group: %s (%d genes)", gid, length(group_list[[gid]])))
    results[[gid]] <- tryCatch(
      score_gene_group(counts_mat, lognorm_mat, det_rate_all,
                       group_list[[gid]], group_id = gid,
                       n_pair_cap = n_pair_cap, seed = seed,
                       timeout_secs = timeout_secs, log_fn = log_fn),
      error = function(e) {
        log_fn(sprintf("  [%s] ERROR: %s", gid, conditionMessage(e)))
        NULL
      }
    )
  }
  results
}

# ── Summarize batch results to data.frame ─────────────────────────────────
summarise_scores <- function(score_list, source_col = "source", source_val = NA_character_) {
  rows <- lapply(score_list, function(r) {
    if (is.null(r)) return(NULL)
    data.frame(
      source            = source_val,
      group_id          = r$group_id,
      n_genes           = r$n_genes,
      n_pairs           = r$n_pairs,
      capped            = r$capped,
      timed_out         = r$timed_out,
      mean_codet_within = r$mean_codet_within,
      mean_codet_null   = r$mean_codet_null,
      gap_codet         = r$gap_codet,
      effect_size_codet = r$effect_size_codet,
      mean_spear_within = r$mean_spear_within,
      mean_spear_null   = r$mean_spear_null,
      gap_spear         = r$gap_spear,
      effect_size_spear = r$effect_size_spear,
      stringsAsFactors  = FALSE
    )
  })
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}
