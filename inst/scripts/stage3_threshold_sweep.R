#!/usr/bin/env Rscript
# Stage 3: edge-threshold selection sweep (FLAG-14)
#
# Evaluates 9 design points using prior-free stability-richness metrics.
# Results are written incrementally (resumable on restart).
#
# Usage:
#   # Full run (Phase 1 + Phase 2b):
#   Rscript inst/scripts/stage3_threshold_sweep.R > logs/stage3_phase2b.log 2>&1
#
#   # Phase 2b only (skip Phase 1 — uses cached edge lists):
#   Rscript inst/scripts/stage3_threshold_sweep.R --phase 2b \
#     > logs/stage3_phase2b.log 2>&1
#
# ENVIRONMENT CONSTRAINT: data.table GForce (frank by=, dt[group,.()],
# unique(by=), dt[bool_expr] on large tables) segfaults on aarch64-darwin.
# All large-table operations use base-R only. fread uses nThread=1L.

suppressPackageStartupMessages({
  library(CoexprArabidopsis)
  library(data.table)
  library(igraph)
})

# Parse --phase flag
.args     <- commandArgs(trailingOnly = TRUE)
.phase_idx <- which(.args == "--phase")
PHASE <- if (length(.phase_idx) && length(.args) > .phase_idx) {
  .args[.phase_idx + 1L]
} else {
  "all"
}

cat(sprintf("Stage 3 threshold sweep  |  phase = %s\n", PHASE))

source("R/coexpr_eval.R")
source("R/stage3_threshold_eval.R")

set.seed(98)

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

SEURAT_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/",
  "Projects/SA_PTI_ETI_single_cell/",
  "SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed/",
  "combined_filtered.rds"
)
DATASET_ID     <- "pathogen_multiome"
STRATUM_VAR    <- "sample2"
SUBCLUSTER_COL <- "sub_clst_rna_20260610"
SYMBOL_MAP_PATH <- "results/pathogen_multiome/symbol_map.csv"

RESULTS_DIR  <- "results/pathogen_multiome/pseudobulk_zscore_spearman"
OUT_DIR      <- "results/pathogen_multiome/stage3_threshold_sweep"
PAIR_SCORES  <- file.path(RESULTS_DIR, "pair_scores_full.csv")

LEVER_A <- c(0.35, 0.40, 0.42, 0.44, 0.46, 0.50)
LEVER_B <- c(30L, 50L, 100L)

BUDGET_SECS  <- 25L * 60L
HARD_TIMEOUT <- 40L * 60L
N_REPS_SPLITHALF <- 5L
N_FOLDS_HELDOUT  <- 5L
N_NULL_PERM      <- 10L
BON3_ID          <- "AT1G08860"

dir.create(OUT_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create("logs",       recursive = TRUE, showWarnings = FALSE)

DENSITY_CSV      <- file.path(OUT_DIR, "density_table.csv")
METRICS_CSV      <- file.path(OUT_DIR, "stage3_metrics.csv")
EDGE_CACHE_RDS   <- file.path(OUT_DIR, "edge_lists_cache.rds")
OBS_CACHE_RDS    <- file.path(OUT_DIR, "obs_normalized_cache.rds")

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

elapsed_secs <- function(t0) as.numeric(proc.time()["elapsed"] - t0["elapsed"])

.safe <- function(expr, fallback, label) {
  tryCatch(expr, error = function(e) {
    message("  [WARN] ", label, ": ", conditionMessage(e))
    fallback
  })
}

append_csv_row <- function(df, path) {
  write_header <- !file.exists(path)
  data.table::fwrite(df, path, append = !write_header, col.names = write_header)
}

stamp <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
  flush.console()
}

# Per-gene rank (1 = highest abs_r) — base-R only (no frank)
.rank_within_gene <- function(abs_r_vec, gene_vec) {
  ord         <- order(gene_vec, -abs_r_vec)
  grp_lengths <- rle(gene_vec[ord])$lengths
  rank_in_grp <- unlist(lapply(grp_lengths, seq_len), use.names = FALSE)
  rank_vec    <- integer(length(abs_r_vec))
  rank_vec[ord] <- rank_in_grp
  rank_vec
}

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 0: adapter check
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("═", 70), "\n", sep = "")
cat("PHASE 0: Adapter compatibility check\n")
cat(strrep("═", 70), "\n")
cat(
  "Stage 3 adapters (obs-point axis — genes×298):\n",
  "  stage3_eval_splithalf       obs-point column split → Jaccard\n",
  "  stage3_eval_effective_rank  SVD of genes×n_pts (masked to visible genes)\n",
  "  stage3_eval_heldout         5-fold CV GBA R² on obs-points\n",
  "  stage3_eval_null_gap        real vs permuted Spearman density\n",
  "  stage3_eval_visible_genes   genes with ≥1 retained edge\n",
  "  stage3_eval_louvain         descriptive only; NOT selection metric\n",
  sep = ""
)
cat("PHASE 0 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 1: Build 9 networks + density table (or load from cache)
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("PHASE 1: Build networks + density table\n")
cat(strrep("═", 70), "\n")

if (PHASE == "2b" && file.exists(EDGE_CACHE_RDS)) {

  stamp("--phase 2b: loading edge lists from cache (skipping CSV read) ...")
  .cache      <- readRDS(EDGE_CACHE_RDS)
  edge_lists  <- .cache$edge_lists
  N_GENES_TOTAL <- .cache$N_GENES_TOTAL
  N_MAX_PAIRS   <- .cache$N_MAX_PAIRS
  density_table <- data.table::fread(DENSITY_CSV)
  rm(.cache)
  stamp(paste("Edge cache loaded:", length(edge_lists), "networks;",
              N_GENES_TOTAL, "genes"))

} else {

  stamp("Reading pair scores (gene_id_A, gene_id_B, z_bar) ...")
  t_p1 <- proc.time()
  dt_pairs <- data.table::fread(PAIR_SCORES,
                                 select   = c("gene_id_A","gene_id_B","z_bar"),
                                 nThread  = 1L)
  stamp(paste("Read", nrow(dt_pairs), "pairs in",
              round(elapsed_secs(t_p1), 1), "s"))

  gA         <- dt_pairs$gene_id_A
  gB         <- dt_pairs$gene_id_B
  mean_abs_r <- tanh(pmin(abs(dt_pairs$z_bar), 9.9))
  rm(dt_pairs); gc(verbose = FALSE)

  all_genes     <- unique(c(gA, gB))
  N_GENES_TOTAL <- length(all_genes)
  N_MAX_PAIRS   <- N_GENES_TOTAL * (N_GENES_TOTAL - 1) / 2
  stamp(paste("Gene universe:", N_GENES_TOTAL, "genes"))

  stamp("Computing per-gene ranks for top-k (base-R) ...")
  rank_as_A <- .rank_within_gene(mean_abs_r, gA)
  rank_as_B <- .rank_within_gene(mean_abs_r, gB)
  stamp("Per-gene ranks done.")

  density_rows <- vector("list", 9L)
  edge_lists   <- vector("list", 9L)

  # Fresh density CSV
  if (file.exists(DENSITY_CSV)) file.remove(DENSITY_CSV)

  .edges_from_idx <- function(idx) {
    data.table(gene_id_A  = gA[idx],
               gene_id_B  = gB[idx],
               mean_abs_r = mean_abs_r[idx])
  }

  cat("\nLever A: global |r| threshold\n")
  for (i in seq_along(LEVER_A)) {
    thr   <- LEVER_A[i]
    idx   <- which(mean_abs_r >= thr)
    edges <- .edges_from_idx(idx)
    ng    <- length(unique(c(edges$gene_id_A, edges$gene_id_B)))
    dens  <- length(idx) / N_MAX_PAIRS
    row   <- data.frame(lever="A_globalr", param=sprintf("%.2f", thr),
                        n_pairs=length(idx), n_genes=ng, density=dens,
                        stringsAsFactors=FALSE)
    density_rows[[i]] <- row;  edge_lists[[i]] <- edges
    append_csv_row(row, DENSITY_CSV)
    cat(sprintf("  |r|>=%.2f  n_pairs=%d  n_genes=%d  density=%.5f\n",
                thr, length(idx), ng, dens))
    flush.console()
  }

  cat("\nLever B: per-gene top-k\n")
  for (j in seq_along(LEVER_B)) {
    k     <- LEVER_B[j]
    idx   <- which(rank_as_A <= k | rank_as_B <= k)
    edges <- .edges_from_idx(idx)
    ng    <- length(unique(c(edges$gene_id_A, edges$gene_id_B)))
    dens  <- length(idx) / N_MAX_PAIRS
    row   <- data.frame(lever="B_topk", param=as.character(k),
                        n_pairs=length(idx), n_genes=ng, density=dens,
                        stringsAsFactors=FALSE)
    density_rows[[6L + j]] <- row;  edge_lists[[6L + j]] <- edges
    append_csv_row(row, DENSITY_CSV)
    cat(sprintf("  top-k=%d   n_pairs=%d  n_genes=%d  density=%.5f\n",
                k, length(idx), ng, dens))
    flush.console()
  }

  density_table <- data.table::rbindlist(density_rows)
  cat("\nDensity table written to", DENSITY_CSV, "\n")

  # Save edge cache for --phase 2b restarts
  saveRDS(list(edge_lists    = edge_lists,
               N_GENES_TOTAL = N_GENES_TOTAL,
               N_MAX_PAIRS   = N_MAX_PAIRS),
          EDGE_CACHE_RDS)
  stamp("Edge cache saved.")

  rm(gA, gB, mean_abs_r, rank_as_A, rank_as_B); gc(verbose = FALSE)
}

cat("\nPHASE 1 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 1b: Load / build normalized obs-point matrix (genes × 298)
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("PHASE 1b: Load obs-point matrix (genes x 298, zscore_gene)\n")
cat(strrep("═", 70), "\n")

if (file.exists(OBS_CACHE_RDS)) {

  stamp("Loading cached obs from", OBS_CACHE_RDS, "...")
  obs <- readRDS(OBS_CACHE_RDS)
  stamp(paste("obs loaded:", nrow(obs$matrix), "genes x",
              ncol(obs$matrix), "obs-points"))

} else {

  if (!file.exists(SEURAT_PATH))
    stop("Seurat object not found: ", SEURAT_PATH)
  if (!file.exists(SYMBOL_MAP_PATH))
    stop("Symbol map not found: ", SYMBOL_MAP_PATH)

  symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)
  stamp("Loading Seurat object (may take 5-10 min) ...")
  t_seurat <- proc.time()

  bundle <- load_seurat(
    seurat_path    = SEURAT_PATH,
    dataset_id     = DATASET_ID,
    stratum_var    = STRATUM_VAR,
    stratum_levels = NULL,          # keep all 13 sample2 levels
    group_var      = SUBCLUSTER_COL,
    symbol_map     = symbol_map,
    min_cells      = 10L
  )
  stamp(paste("Bundle loaded in", round(elapsed_secs(t_seurat) / 60, 1), "min:",
              nrow(bundle$counts), "genes x", ncol(bundle$counts), "cells"))

  stamp("Building obs_subcluster ...")
  obs_raw    <- obs_subcluster(bundle, group_col = SUBCLUSTER_COL)
  obs        <- obs_raw
  obs$matrix <- normalize_obs(obs_raw, method = "zscore_gene")
  rm(bundle, obs_raw); gc(verbose = FALSE)

  stamp(paste("obs:", nrow(obs$matrix), "genes x",
              ncol(obs$matrix), "obs-points (zscore_gene)"))

  # Cache to avoid re-loading Seurat on restart
  saveRDS(obs, OBS_CACHE_RDS)
  stamp(paste("obs cached to", OBS_CACHE_RDS))
}

cat("\nPHASE 1b complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 2: Prior-free evaluation at each of the 9 design points
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("PHASE 2: Prior-free evaluation (obs-point axis)\n")
cat(strrep("═", 70), "\n")

# Assemble design points
design_points <- vector("list", 9L)
for (i in seq_along(LEVER_A))
  design_points[[i]] <- list(lever       = "A_globalr",
                              param       = sprintf("%.2f", LEVER_A[i]),
                              threshold_r = LEVER_A[i],
                              top_k       = NULL,
                              edges       = edge_lists[[i]])
for (j in seq_along(LEVER_B))
  design_points[[6L + j]] <- list(lever       = "B_topk",
                                   param       = as.character(LEVER_B[j]),
                                   threshold_r = NULL,
                                   top_k       = LEVER_B[j],
                                   edges       = edge_lists[[6L + j]])

# Resumability: skip points already in METRICS_CSV
already_done <- character(0)
if (file.exists(METRICS_CSV)) {
  prev <- data.table::fread(METRICS_CSV)
  if (nrow(prev) > 0L && all(c("lever","param") %in% names(prev)))
    already_done <- paste0(prev$lever, "_", prev$param)
  stamp(paste("Resuming — skipping", length(already_done), "completed point(s)"))
}

for (pt in design_points) {
  pt_key <- paste0(pt$lever, "_", pt$param)

  if (pt_key %in% already_done) {
    cat("  SKIP", pt_key, "(already in CSV)\n"); flush.console()
    next
  }

  cat("\n", strrep("-", 60), "\n", sep = "")
  stamp(paste("POINT", pt_key, "| n_pairs =", nrow(pt$edges)))
  t_pt     <- proc.time()
  notes    <- character(0)
  timedout <- FALSE   # single flag; avoids exists() checks across loop body

  # NA defaults — overwritten if metrics succeed
  vis_res  <- data.frame(n_visible=NA_integer_, n_total=N_GENES_TOTAL, frac_visible=NA_real_)
  eff_res  <- data.frame(eff_rank=NA_real_, n_visible=NA_integer_, n_points=NA_integer_)
  null_res <- data.frame(null_gap=NA_real_, real_frac=NA_real_, perm_frac_mean=NA_real_)
  pred_res <- data.frame(heldout_r2=NA_real_)
  sh_res   <- data.frame(splithalf_jaccard=NA_real_, splithalf_pearson=NA_real_,
                          splithalf_jaccard_sd=NA_real_, n_reps=0L)
  lou_res  <- data.frame(n_modules=NA_integer_, grey_rate=NA_real_,
                          median_module_size=NA_real_)
  n_reps_sh <- N_REPS_SPLITHALF

  # ── 1. Visible genes (always run) ───────────────────────────────────────────
  stamp("  computing visible_genes ...")
  vis_res <- .safe(
    stage3_eval_visible_genes(pt$edges, n_total = N_GENES_TOTAL),
    vis_res, "visible_genes"
  )
  cat("    n_visible =", vis_res$n_visible, "\n"); flush.console()

  # ── 2. Effective rank ────────────────────────────────────────────────────────
  if (!timedout) {
    stamp("  computing eff_rank ...")
    eff_res <- .safe(
      stage3_eval_effective_rank(obs, pt$edges),
      eff_res, "eff_rank"
    )
    cat("    eff_rank =", round(eff_res$eff_rank, 3), "\n"); flush.console()
    if (elapsed_secs(t_pt) > HARD_TIMEOUT) {
      notes <- c(notes, paste0("TIMED OUT after eff_rank (",
                                round(elapsed_secs(t_pt)/60,1), " min)"))
      timedout <- TRUE
      stamp("  HARD TIMEOUT — skipping null_gap + heldout + splithalf")
    }
  }

  # ── 3. Null gap ──────────────────────────────────────────────────────────────
  if (!timedout) {
    if (elapsed_secs(t_pt) > BUDGET_SECS) {
      notes <- c(notes, "null_gap skipped: >25 min budget")
      n_reps_sh <- 3L
      notes <- c(notes, "splithalf reduced to n_reps=3")
      stamp("  SKIPPING null_gap (over 25-min budget); splithalf→3 reps")
    } else {
      stamp(paste("  computing null_gap (", N_NULL_PERM, "perms) ..."))
      null_res <- .safe(
        stage3_eval_null_gap(obs, pt$edges,
                             min_abs_r = pt$threshold_r, top_k = pt$top_k,
                             n_perm = N_NULL_PERM, seed = 98L),
        null_res, "null_gap"
      )
      cat("    null_gap =", round(null_res$null_gap, 3), "\n"); flush.console()
    }
    if (elapsed_secs(t_pt) > HARD_TIMEOUT) {
      notes <- c(notes, paste0("TIMED OUT after null_gap (",
                                round(elapsed_secs(t_pt)/60,1), " min)"))
      timedout <- TRUE
      stamp("  HARD TIMEOUT — skipping heldout + splithalf")
    }
  }

  # ── 4. Heldout predictivity ──────────────────────────────────────────────────
  if (!timedout) {
    stamp("  computing heldout_r2 ...")
    pred_res <- .safe(
      stage3_eval_heldout(obs, pt$edges,
                          min_abs_r = pt$threshold_r, top_k = pt$top_k,
                          n_folds = N_FOLDS_HELDOUT, seed = 98L),
      pred_res, "heldout"
    )
    cat("    heldout_r2 =", round(pred_res$heldout_r2, 4), "\n"); flush.console()
    if (elapsed_secs(t_pt) > HARD_TIMEOUT) {
      notes <- c(notes, paste0("TIMED OUT after heldout (",
                                round(elapsed_secs(t_pt)/60,1), " min)"))
      timedout <- TRUE
      stamp("  HARD TIMEOUT — skipping splithalf")
    }
  }

  # ── 5. Split-half ────────────────────────────────────────────────────────────
  if (!timedout) {
    stamp(paste("  computing split-half (n_reps =", n_reps_sh, ") ..."))
    sh_res <- .safe(
      stage3_eval_splithalf(obs, pt$edges,
                            min_abs_r = pt$threshold_r, top_k = pt$top_k,
                            n_reps = n_reps_sh, seed = 98L),
      sh_res, "splithalf"
    )
    cat("    splithalf_jaccard =", round(sh_res$splithalf_jaccard, 4),
        "(pearson =", round(sh_res$splithalf_pearson, 4), ")\n")
    flush.console()
    if (elapsed_secs(t_pt) > HARD_TIMEOUT) {
      notes <- c(notes, paste0("TIMED OUT after splithalf (",
                                round(elapsed_secs(t_pt)/60,1), " min)"))
      timedout <- TRUE
      stamp("  HARD TIMEOUT — skipping Louvain")
    }
  }

  # ── 6. Louvain (descriptive) ─────────────────────────────────────────────────
  if (!timedout) {
    stamp("  computing Louvain modules (descriptive) ...")
    lou_res <- .safe(
      stage3_eval_louvain(pt$edges, seed = 98L),
      lou_res, "louvain"
    )
    cat("    n_modules =", lou_res$n_modules,
        "| grey_rate =", round(lou_res$grey_rate, 3), "\n")
    flush.console()
  }

  eval_secs <- round(elapsed_secs(t_pt), 1)
  cat("  Total elapsed:", round(eval_secs/60, 2), "min\n")

  # ── Write result row ──────────────────────────────────────────────────────────
  den_row <- density_table[density_table$lever == pt$lever &
                             density_table$param == pt$param, ]
  out_row <- data.frame(
    lever                      = pt$lever,
    param                      = pt$param,
    density                    = if (nrow(den_row)) den_row$density[1L] else NA_real_,
    n_genes                    = vis_res$n_visible,
    splithalf_jaccard          = sh_res$splithalf_jaccard,
    splithalf_pearson          = sh_res$splithalf_pearson,
    eff_rank                   = eff_res$eff_rank,
    heldout_r2                 = pred_res$heldout_r2,
    null_gap                   = null_res$null_gap,
    visible_genes              = vis_res$n_visible,
    louvain_n_modules          = lou_res$n_modules,
    louvain_grey_rate          = lou_res$grey_rate,
    louvain_median_module_size = lou_res$median_module_size,
    eval_seconds               = eval_secs,
    notes                      = paste(notes, collapse = "; "),
    stringsAsFactors           = FALSE
  )
  append_csv_row(out_row, METRICS_CSV)
  stamp(paste("POINT", pt_key, "complete — jaccard =",
              round(sh_res$splithalf_jaccard, 4),
              "eff_rank =", round(eff_res$eff_rank, 3)))

  rm(vis_res, eff_res, null_res, pred_res, sh_res, lou_res, out_row, den_row,
     timedout, n_reps_sh, notes, eval_secs)
  gc(verbose = FALSE)
}

cat("\nPHASE 2 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3: Pareto analysis + recommendation
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("PHASE 3: Pareto analysis + recommendation\n")
cat(strrep("═", 70), "\n")

metrics <- data.table::fread(METRICS_CSV)

cat("\nMetrics table:\n")
print(metrics[, .(lever, param, density, splithalf_jaccard, eff_rank,
                   heldout_r2, visible_genes, louvain_n_modules)])
flush.console()

# Pareto front: splithalf_jaccard (stability) vs eff_rank + heldout_r2 (richness)
is_dominated <- function(dt) {
  sh <- dt$splithalf_jaccard
  er <- dt$eff_rank
  hr <- dt$heldout_r2
  n  <- nrow(dt)
  dom <- logical(n)
  for (i in seq_len(n)) {
    if (is.na(sh[i]) || is.na(er[i])) { dom[i] <- TRUE; next }
    hr_i <- if (is.na(hr[i])) -Inf else hr[i]
    for (j in seq_len(n)) {
      if (i == j || is.na(sh[j]) || is.na(er[j])) next
      hr_j <- if (is.na(hr[j])) -Inf else hr[j]
      if (sh[j] >= sh[i] && er[j] >= er[i] && hr_j >= hr_i &&
          (sh[j] > sh[i] || er[j] > er[i] || hr_j > hr_i)) {
        dom[i] <- TRUE; break
      }
    }
  }
  dom
}

metrics[, dominated := is_dominated(.SD)]
pareto_pts <- metrics[dominated == FALSE]

cat("\nPareto-front design points (splithalf × eff_rank × heldout_r2):\n")
print(pareto_pts[, .(lever, param, density, splithalf_jaccard, eff_rank, heldout_r2)])
flush.console()

# Lever A vs Lever B at matched density
cat("\nLever A vs Lever B at matched density:\n")
dens_A <- metrics[lever == "A_globalr"]
dens_B <- metrics[lever == "B_topk"]
lever_compare <- list()
for (bi in seq_len(nrow(dens_B))) {
  d_B   <- dens_B$density[bi]
  diffs <- abs(dens_A$density - d_B)
  best  <- which.min(diffs)
  if (!length(best)) next
  sh_A <- dens_A$splithalf_jaccard[best]; sh_B <- dens_B$splithalf_jaccard[bi]
  er_A <- dens_A$eff_rank[best];         er_B <- dens_B$eff_rank[bi]
  hr_A <- dens_A$heldout_r2[best];       hr_B <- dens_B$heldout_r2[bi]

  v_sh <- function(a, b) if (!is.na(a) && !is.na(b)) ifelse(b>a,"B better",ifelse(b<a,"A better","tie")) else "NA"
  cat(sprintf("  top-k=%s (%.5f) vs |r|>=%.2f (%.5f):\n",
              dens_B$param[bi], d_B, as.numeric(dens_A$param[best]),
              dens_A$density[best]))
  cat(sprintf("    splithalf: B=%.4f  A=%.4f  → %s\n", sh_B, sh_A, v_sh(sh_A, sh_B)))
  cat(sprintf("    eff_rank:  B=%.3f  A=%.3f  → %s\n", er_B, er_A, v_sh(er_A, er_B)))
  cat(sprintf("    heldout:   B=%.4f  A=%.4f  → %s\n",
              ifelse(is.na(hr_B),-99,hr_B), ifelse(is.na(hr_A),-99,hr_A),
              v_sh(hr_A, hr_B)))
  lever_compare[[bi]] <- data.frame(k_B=dens_B$param[bi], d_B=d_B,
    thr_A=dens_A$param[best], d_A=dens_A$density[best],
    sh_A=sh_A, sh_B=sh_B, er_A=er_A, er_B=er_B, hr_A=hr_A, hr_B=hr_B,
    v_sh=v_sh(sh_A,sh_B), v_er=v_sh(er_A,er_B), v_hr=v_sh(hr_A,hr_B),
    stringsAsFactors=FALSE)
  flush.console()
}

# Recommendation
if (nrow(pareto_pts) > 0L) {
  pp2 <- pareto_pts[!is.na(splithalf_jaccard) & !is.na(eff_rank)]
  if (nrow(pp2) == 0L) pp2 <- pareto_pts
  hr_rng <- range(pp2$heldout_r2, na.rm = TRUE)
  er_rng <- range(pp2$eff_rank,   na.rm = TRUE)
  pp2[, composite := splithalf_jaccard +
        0.3 * (if (diff(er_rng) > 0) (eff_rank - er_rng[1]) / diff(er_rng) else 0) +
        0.2 * (if (diff(hr_rng) > 0) (ifelse(is.na(heldout_r2), hr_rng[1], heldout_r2) -
                                       hr_rng[1]) / diff(hr_rng) else 0)]
  rec <- pp2[which.max(composite)]
} else {
  rec <- metrics[which.max(splithalf_jaccard +
                            0.3 * ifelse(is.na(eff_rank), 0, eff_rank) +
                            0.2 * ifelse(is.na(heldout_r2), 0, heldout_r2))]
}

cat("\n--- RECOMMENDATION (prior-free metrics only) ---\n")
cat("Recommended: lever =", rec$lever, "  param =", rec$param, "\n")
cat("  density          =", round(rec$density, 5), "\n")
cat("  splithalf_jaccard=", round(rec$splithalf_jaccard, 4), " (stability)\n")
cat("  eff_rank         =", round(rec$eff_rank, 3), " (richness)\n")
cat("  heldout_r2       =", round(rec$heldout_r2, 4), " (richness)\n")
cat("  null_gap         =", round(rec$null_gap, 3), "\n")
cat("  visible_genes    =", rec$visible_genes, "\n")
cat("NOTE: recommendation from prior-free metrics only.",
    "Final threshold call is the user's.\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3b: Post-hoc sanity at recommended point (NOT selection)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n--- POST-HOC SANITY (recommended point; NOT a selection input) ---\n")

rec_idx   <- which(sapply(design_points,
                           function(p) paste0(p$lever,"_",p$param)) ==
                     paste0(rec$lever,"_",rec$param))
edges_rec <- if (length(rec_idx)) {
  design_points[[rec_idx[1L]]]$edges
} else {
  data.table(gene_id_A=character(), gene_id_B=character(), mean_abs_r=numeric())
}

net_genes <- unique(c(edges_rec$gene_id_A, edges_rec$gene_id_B))

# BON3
bon3_partners <- unique(c(
  edges_rec$gene_id_B[edges_rec$gene_id_A == BON3_ID],
  edges_rec$gene_id_A[edges_rec$gene_id_B == BON3_ID]
))
cat("BON3 (AT1G08860):\n")
cat("  in network:", BON3_ID %in% net_genes, "\n")
cat("  n partners:", length(bon3_partners), "\n")

# WRKY
sym_paths <- c("results/pathogen_multiome/symbol_map.csv", "results/symbol_map.csv")
wrky_ids  <- character(0)
for (sp in sym_paths) {
  if (file.exists(sp)) {
    sym  <- data.table::fread(sp)
    wrky_ids <- sym[[1L]][grepl("^WRKY", sym[[2L]])]
    break
  }
}
n_wrky_net <- sum(wrky_ids %in% net_genes)
cat("WRKY family (n=", length(wrky_ids), "):\n", sep = "")
cat("  in network:", n_wrky_net, "/", length(wrky_ids),
    sprintf("(%.1f%%)\n", 100 * n_wrky_net / max(length(wrky_ids), 1L)))
cat("(Post-hoc sanity only — NOT used in threshold selection)\n")

# Louvain kME at recommended point
cat("\n--- kME POST-HOC (recommended point; NOT selection) ---\n")
kme_df     <- NULL
b3_kme_row <- NULL
wrky_kme_df <- NULL

tryCatch({
  ns_rec <- .s3_net_subset(obs, edges_rec)

  # Build Spearman fingerprint matrix for the recommended point (full obs)
  mat_fp   <- obs$matrix[rownames(obs$matrix) %in% net_genes, , drop = FALSE]
  gids     <- rownames(mat_fp)

  # Louvain for module membership
  g_rec  <- igraph::graph_from_data_frame(
    data.frame(gene_id_A = edges_rec$gene_id_A,
               gene_id_B = edges_rec$gene_id_B,
               weight    = edges_rec$mean_abs_r,
               stringsAsFactors = FALSE),
    directed = FALSE
  )
  set.seed(98L)
  cl_rec   <- igraph::cluster_louvain(g_rec, weights = igraph::E(g_rec)$weight)
  memb_vec <- igraph::membership(cl_rec)

  kme_df <- data.frame(gene_id = gids,
                        module  = as.integer(memb_vec[gids]),
                        kme     = NA_real_,
                        stringsAsFactors = FALSE)
  gene_idx <- setNames(seq_len(nrow(kme_df)), kme_df$gene_id)

  # kME = correlation of gene fingerprint (mean neighbor r per obs-point)
  # with module mean fingerprint. With 298 obs-points this is meaningful.
  # Fingerprint: for each gene g, its fingerprint is its expression profile
  # across 298 obs-points (from obs$matrix). Module mean = colMeans of members.
  for (mod_id in sort(unique(kme_df$module))) {
    gi <- kme_df$gene_id[kme_df$module == mod_id]
    if (length(gi) < 2L) next
    mod_mat  <- mat_fp[gi, , drop = FALSE]
    mod_mean <- colMeans(mod_mat)
    if (stats::sd(mod_mean) < 1e-10) next
    for (g in gi)
      kme_df$kme[gene_idx[g]] <- stats::cor(mat_fp[g, ], mod_mean)
  }

  b3_kme_row  <- kme_df[kme_df$gene_id == BON3_ID, ]
  wrky_kme_df <- kme_df[kme_df$gene_id %in% wrky_ids, ]

  cat("BON3 kME:\n")
  if (nrow(b3_kme_row) == 1L && !is.na(b3_kme_row$kme))
    cat(sprintf("  module=%d  kME=%.4f\n", b3_kme_row$module, b3_kme_row$kme))
  else
    cat("  BON3 not in network or singleton module (kME=NA)\n")

  n_w_net  <- nrow(wrky_kme_df)
  n_w_asgn <- sum(!is.na(wrky_kme_df$kme))
  cat(sprintf("WRKY: %d in network | %d non-singleton (kME assigned)\n",
              n_w_net, n_w_asgn))
  if (n_w_asgn > 0L) {
    kvals <- sort(wrky_kme_df$kme[!is.na(wrky_kme_df$kme)], decreasing = TRUE)
    cat(sprintf("  kME: min=%.3f Q1=%.3f median=%.3f Q3=%.3f max=%.3f\n",
                min(kvals), stats::quantile(kvals,.25), stats::median(kvals),
                stats::quantile(kvals,.75), max(kvals)))
    top5 <- wrky_kme_df[order(-wrky_kme_df$kme), ][seq_len(min(5L,n_w_asgn)),]
    cat("  Top-5 WRKY by kME:\n")
    for (i in seq_len(nrow(top5)))
      cat(sprintf("    %s  mod=%d  kME=%.4f\n",
                  top5$gene_id[i], top5$module[i], top5$kme[i]))
  }
}, error = function(e) {
  message("  [WARN] kME computation failed: ", conditionMessage(e))
})
cat("(kME post-hoc only — NOT a selection input)\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3c: Write STAGE3_FINDINGS.md
# ─────────────────────────────────────────────────────────────────────────────

bon3_kme_str <- if (!is.null(b3_kme_row) && nrow(b3_kme_row)==1 && !is.na(b3_kme_row$kme)) {
  sprintf("module=%d  kME=%.4f", b3_kme_row$module, b3_kme_row$kme)
} else {
  "kME=NA (not in network or singleton module)"
}

wrky_kme_str <- if (!is.null(wrky_kme_df) && sum(!is.na(wrky_kme_df$kme)) > 0L) {
  kv   <- sort(wrky_kme_df$kme[!is.na(wrky_kme_df$kme)], decreasing=TRUE)
  top5 <- wrky_kme_df[order(-wrky_kme_df$kme),][seq_len(min(5L,length(kv))),]
  t5s  <- paste(sprintf("%s(mod=%d,kME=%.3f)", top5$gene_id, top5$module, top5$kme),
                collapse="; ")
  sprintf("n_in_net=%d assigned=%d; kME median=%.3f [%.3f,%.3f]; top-5: %s",
          nrow(wrky_kme_df), length(kv), median(kv), min(kv), max(kv), t5s)
} else "kME unavailable"

.fmt_md_table <- function(dt) {
  cols <- names(dt)
  hdr  <- paste("|", paste(cols, collapse=" | "), "|")
  sep  <- paste("|", paste(rep("---", length(cols)), collapse=" | "), "|")
  rows <- apply(as.data.frame(dt), 1L,
                function(r) paste("|", paste(r, collapse=" | "), "|"))
  paste(c(hdr, sep, rows), collapse="\n")
}

pareto_tbl <- if (nrow(pareto_pts) > 0L) {
  .fmt_md_table(pareto_pts[, .(lever, param,
                                density        = round(density,5),
                                splithalf      = round(splithalf_jaccard,4),
                                eff_rank       = round(eff_rank,3),
                                heldout_r2     = round(heldout_r2,4))])
} else "No non-dominated Pareto points (check for NA metrics)."

lever_md <- paste(sapply(lever_compare, function(lc) {
  if (is.null(lc)) return("")
  paste0("**top-k=", lc$k_B, " vs |r|>=", lc$thr_A, "**\n",
         "- splithalf: B=", round(lc$sh_B,4), "  A=", round(lc$sh_A,4),
         "  → ", lc$v_sh, "\n",
         "- eff_rank:  B=", round(lc$er_B,3), "  A=", round(lc$er_A,3),
         "  → ", lc$v_er, "\n",
         "- heldout:   B=", round(ifelse(is.na(lc$hr_B),-99,lc$hr_B),4),
         "  A=", round(ifelse(is.na(lc$hr_A),-99,lc$hr_A),4),
         "  → ", lc$v_hr)
}), collapse="\n\n")

findings_md <- c(
  "# Stage 3: Edge-Threshold Selection Findings (Phase 2b — obs-point basis)",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Metrics use the obs-point axis (genes × 298 pseudobulk profiles), NOT the",
  "4-condition fingerprint from the prior invalid Phase 2.",
  "",
  "## 1. Density Table (9 design points)",
  "", .fmt_md_table(density_table), "",
  "## 2. Prior-Free Metrics",
  "", .fmt_md_table(metrics[, .(lever, param,
                                  density=round(density,5),
                                  splithalf=round(splithalf_jaccard,4),
                                  eff_rank=round(eff_rank,3),
                                  null_gap=round(null_gap,2),
                                  heldout_r2=round(heldout_r2,4),
                                  visible_genes,
                                  louvain_n_modules,
                                  louvain_grey_rate=round(louvain_grey_rate,3),
                                  eval_seconds, notes)]), "",
  "## 3. Pareto Front (splithalf × eff_rank × heldout_r2)", "",
  pareto_tbl, "",
  "## 4. Recommended Design Point",
  "",
  paste0("**", rec$lever, "  param = ", rec$param, "**"), "",
  "| metric | value |", "| --- | --- |",
  paste0("| density | ", round(rec$density,5), " |"),
  paste0("| splithalf_jaccard | ", round(rec$splithalf_jaccard,4), " |"),
  paste0("| eff_rank | ", round(rec$eff_rank,3), " |"),
  paste0("| heldout_r2 | ", round(rec$heldout_r2,4), " |"),
  paste0("| null_gap | ", round(rec$null_gap,3), " |"),
  paste0("| visible_genes | ", rec$visible_genes, " |"),
  "",
  "**Rationale:** Pareto-dominant on stability (splithalf_jaccard) and richness",
  "(eff_rank, heldout_r2) using the 298-obs-point prior-free harness.",
  "No gold-standard gene sets used in selection.",
  "",
  "> **Final threshold choice is the user's.**",
  "",
  "## 5. Lever A vs Lever B at Matched Density", "",
  lever_md, "",
  "## 6. Post-Hoc Sanity at Recommended Point", "",
  "*(NOT used for selection.)*", "",
  paste0("**BON3 (AT1G08860):** in network=", BON3_ID %in% net_genes,
         ";  n partners=", length(bon3_partners), ";  ", bon3_kme_str),
  "",
  paste0("**WRKY (n=", length(wrky_ids), "):** ", n_wrky_net, " in network (",
         round(100*n_wrky_net/max(length(wrky_ids),1),1), "%);  ", wrky_kme_str),
  ""
)

writeLines(findings_md, file.path(OUT_DIR, "STAGE3_FINDINGS.md"))
cat("\nFindings written to", file.path(OUT_DIR, "STAGE3_FINDINGS.md"), "\n")
cat("PHASE 3 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# FINAL REPORT
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("FINAL REPORT\n")
cat(strrep("═", 70), "\n\n")

cat("1. DENSITY TABLE (9 rows — unchanged from Phase 1)\n")
print(density_table)

cat("\n2. METRICS TABLE\n")
print(metrics[, .(lever, param, density, splithalf_jaccard, eff_rank,
                   heldout_r2, null_gap, visible_genes, eval_seconds)])

cat("\n3. PARETO FRONT + RECOMMENDATION\n")
if (nrow(pareto_pts) > 0L) {
  print(pareto_pts[, .(lever, param, density, splithalf_jaccard, eff_rank, heldout_r2)])
} else {
  cat("  (no non-dominated points — check for NA metrics)\n")
}
cat("Recommended:", rec$lever, " param =", rec$param, "\n")
cat("  stability (splithalf_jaccard):", round(rec$splithalf_jaccard, 4), "\n")
cat("  richness  (eff_rank):",          round(rec$eff_rank, 3), "\n")
cat("  richness  (heldout_r2):",        round(rec$heldout_r2, 4), "\n")

cat("\n4. LEVER A vs LEVER B AT MATCHED DENSITY\n")
for (lc in lever_compare) {
  if (is.null(lc)) next
  cat(sprintf("  top-k=%s vs |r|>=%.2f: splithalf %s, eff_rank %s, heldout %s\n",
              lc$k_B, as.numeric(lc$thr_A), lc$v_sh, lc$v_er, lc$v_hr))
}

cat("\n5. BON3 + WRKY POST-HOC SANITY (NOT selection criteria)\n")
cat("   BON3 in network:", BON3_ID %in% net_genes,
    "| n partners:", length(bon3_partners), "\n")
cat("   BON3 kME:", bon3_kme_str, "\n")
cat("   WRKY in network:", n_wrky_net, "/", length(wrky_ids), "\n")
cat("   WRKY kME:", wrky_kme_str, "\n")

cat("\nStage 3 Phase 2b sweep complete.\n")
flush.console()
