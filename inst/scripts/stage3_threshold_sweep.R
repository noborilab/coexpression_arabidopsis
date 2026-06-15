#!/usr/bin/env Rscript
# Stage 3: edge-threshold selection benchmark (FLAG-14)
#
# Evaluates 9 design points (6 global-|r| + 3 per-gene top-k) using
# prior-free stability-richness metrics. Writes results incrementally so
# the run is resumable on restart: on re-launch, rows already present in
# stage3_metrics.csv are skipped.
#
# Usage:
#   cd <repo_root>
#   nohup Rscript inst/scripts/stage3_threshold_sweep.R > stage3_sweep.log 2>&1 &

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})
source("R/coexpr_eval.R")
source("R/stage3_threshold_eval.R")

set.seed(98)

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

RESULTS_DIR  <- "results/pathogen_multiome/pseudobulk_zscore_spearman"
OUT_DIR      <- "results/pathogen_multiome/stage3_threshold_sweep"
PAIR_SCORES  <- file.path(RESULTS_DIR, "pair_scores_full.csv")

COR_FILES <- setNames(
  file.path(RESULTS_DIR,
            paste0("cor_", c("Mock","DC3000","AvrRpt2","AvrRpm1"), ".rds")),
  c("r_Mock","r_DC3000","r_AvrRpt2","r_AvrRpm1")
)
COR_COLS <- names(COR_FILES)

N_GENES_TOTAL <- 11010L
N_MAX_PAIRS   <- N_GENES_TOTAL * (N_GENES_TOTAL - 1L) / 2L   # 60 604 545

LEVER_A <- c(0.35, 0.40, 0.42, 0.44, 0.46, 0.50)
LEVER_B <- c(30L, 50L, 100L)

BUDGET_SECS  <- 25L * 60L   # 25 min per point
HARD_TIMEOUT <- 40L * 60L   # 40 min absolute ceiling per point
N_NULL_PERM  <- 10L
BON3_ID      <- "AT1G08860"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DENSITY_CSV <- file.path(OUT_DIR, "density_table.csv")
METRICS_CSV <- file.path(OUT_DIR, "stage3_metrics.csv")

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

stamp <- function(msg) {
  cat(format(Sys.time(), "[%H:%M:%S]"), msg, "\n")
  flush.console()
}

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 0: adapter compatibility check
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("═", 70), "\n", sep = "")
cat("PHASE 0: Adapter compatibility check\n")
cat(strrep("═", 70), "\n")

cat(
  "coexpr_eval.R metrics and Stage 3 status:\n",
  "  eval_splithalf              ADAPTED  stage3_eval_splithalf\n",
  "                              -> 3 unique 2+2 condition splits\n",
  "  eval_effective_rank         ADAPTED  stage3_eval_effective_rank\n",
  "                              -> SVD of genes×4 fingerprint (rank ≤ 4)\n",
  "  eval_null_gap               ADAPTED  stage3_eval_null_gap\n",
  "                              -> fingerprint corr. threshold=0.3, max 2000 genes\n",
  "  eval_heldout_predictivity   ADAPTED  stage3_eval_heldout_predictivity\n",
  "                              -> leave-one-condition-out GBA R²\n",
  "  eval_visible_genes          ADAPTED  stage3_eval_visible_genes\n",
  "                              -> genes with ≥1 retained edge\n",
  "  stage3_eval_louvain         NEW (descriptive only; NOT selection metric)\n",
  "  eval_downsample_depth       SKIPPED  (requires raw count matrix)\n",
  "  eval_downsample_cells       SKIPPED  (requires raw cell matrix)\n",
  "  eval_depth_leakage          SKIPPED  (depth pre-adjusted)\n",
  sep = ""
)

cat("PHASE 0 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 1: Build 9 networks + record density
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("PHASE 1: Build networks + density table\n")
cat(strrep("═", 70), "\n")

stamp("Reading pair scores (gene_id_A, gene_id_B, z_bar) ...")
t_p1 <- proc.time()
dt_pairs <- fread(PAIR_SCORES, select = c("gene_id_A","gene_id_B","z_bar"))
stamp(paste("Read", nrow(dt_pairs), "pairs in", round(elapsed_secs(t_p1),1), "s"))

# mean |r| = tanh(|z_bar|); cap z_bar at 9.9 to avoid Inf from R_score=1
dt_pairs[, abs_r := tanh(pmin(abs(z_bar), 9.9))]

# Per-gene ranks (once, reused for all Lever B design points)
stamp("Computing per-gene ranks for top-k ...")
dt_pairs[, rank_as_A := frank(-abs_r, ties.method = "first"), by = gene_id_A]
dt_pairs[, rank_as_B := frank(-abs_r, ties.method = "first"), by = gene_id_B]
stamp("Per-gene ranks done.")

density_rows <- vector("list", 9L)
edge_lists   <- vector("list", 9L)

# ── Lever A ──────────────────────────────────────────────────────────────────
cat("\nLever A: global |r| threshold\n")
for (i in seq_along(LEVER_A)) {
  thr   <- LEVER_A[i]
  edges <- dt_pairs[abs_r >= thr, .(gene_id_A, gene_id_B, abs_r)]
  npairs <- nrow(edges)
  ngenes <- length(unique(c(edges$gene_id_A, edges$gene_id_B)))
  dens   <- npairs / N_MAX_PAIRS

  density_rows[[i]] <- data.frame(lever="A_globalr", param=sprintf("%.2f",thr),
                                   n_pairs=npairs, n_genes=ngenes, density=dens,
                                   stringsAsFactors=FALSE)
  edge_lists[[i]] <- edges
  cat(sprintf("  |r|>=%.2f  n_pairs=%d  n_genes=%d  density=%.5f\n",
              thr, npairs, ngenes, dens))
  flush.console()
}

# ── Lever B ──────────────────────────────────────────────────────────────────
cat("\nLever B: per-gene top-k (union)\n")
for (j in seq_along(LEVER_B)) {
  k     <- LEVER_B[j]
  edges <- dt_pairs[rank_as_A <= k | rank_as_B <= k,
                    .(gene_id_A, gene_id_B, abs_r)]
  npairs <- nrow(edges)
  ngenes <- length(unique(c(edges$gene_id_A, edges$gene_id_B)))
  dens   <- npairs / N_MAX_PAIRS

  density_rows[[6L + j]] <- data.frame(lever="B_topk", param=as.character(k),
                                        n_pairs=npairs, n_genes=ngenes, density=dens,
                                        stringsAsFactors=FALSE)
  edge_lists[[6L + j]] <- edges
  cat(sprintf("  top-k=%d   n_pairs=%d  n_genes=%d  density=%.5f\n",
              k, npairs, ngenes, dens))
  flush.console()
}

density_table <- rbindlist(density_rows)
fwrite(density_table, DENSITY_CSV)
cat("\nDensity table written to", DENSITY_CSV, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# Pre-extract per-condition r for candidate pairs
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("-", 60), "\n", sep = "")
stamp("Building candidate pair tables for per-condition extraction ...")

# Global-r candidates: abs_r >= 0.27 (buffer 0.08 below minimum threshold 0.35)
cand_global <- dt_pairs[abs_r >= 0.27, .(gene_id_A, gene_id_B, abs_r)]
stamp(paste("  Global-r candidates (abs_r>=0.27):", nrow(cand_global)))

# Top-k candidates: top-200 per gene (covers top-100 with buffer for splits)
cand_topk <- dt_pairs[rank_as_A <= 200L | rank_as_B <= 200L,
                       .(gene_id_A, gene_id_B, abs_r)]
stamp(paste("  Top-200 candidates per gene:", nrow(cand_topk)))

# Union of both (removes duplicates)
cand_all <- unique(
  rbind(cand_global, cand_topk),
  by = c("gene_id_A","gene_id_B")
)
stamp(paste("  Union candidates:", nrow(cand_all)))

# Free large tables before loading cor matrices
rm(dt_pairs, cand_global)
gc(verbose = FALSE)
stamp("Freed dt_pairs. Starting per-condition r extraction ...")

# Gene name → integer index (from first cor matrix; all share the same rownames)
gene_names_universe <- NULL

for (cond_col in COR_COLS) {
  cond_name <- sub("^r_", "", cond_col)
  stamp(paste("  Loading cor_", cond_name, ".rds ...", sep=""))
  t_load  <- proc.time()
  cor_mat <- readRDS(COR_FILES[cond_col])

  if (is.null(gene_names_universe))
    gene_names_universe <- rownames(cor_mat)

  gidx <- setNames(seq_along(gene_names_universe), gene_names_universe)
  i_A  <- gidx[cand_all$gene_id_A]
  i_B  <- gidx[cand_all$gene_id_B]

  cand_all[[cond_col]] <- cor_mat[cbind(i_A, i_B)]

  rm(cor_mat); gc(verbose = FALSE)
  stamp(paste("  Done in", round(elapsed_secs(t_load),1), "s"))
}

# Derive top-200 subset from cand_all for Lever B split-half
# (avoids loading cand_topk separately with per-cond r; re-rank within cand_all)
stamp("Building top-200-per-gene subset from cand_all for Lever B split-half ...")
cand_all[, .rank_A := frank(-abs_r, ties.method = "first"), by = gene_id_A]
cand_all[, .rank_B := frank(-abs_r, ties.method = "first"), by = gene_id_B]
cand_topk_sh <- cand_all[.rank_A <= 200L | .rank_B <= 200L]
cand_all[, c(".rank_A", ".rank_B") := NULL]
stamp(paste("  Top-200 split-half candidates:", nrow(cand_topk_sh)))

# Persist to disk (for debugging / resume)
saveRDS(cand_all,      file.path(OUT_DIR, "cand_extended.rds"),    compress = FALSE)
saveRDS(cand_topk_sh,  file.path(OUT_DIR, "cand_topk_sh.rds"),     compress = FALSE)
rm(cand_topk); gc(verbose = FALSE)
stamp("Pre-extraction complete.")

cat("\nPHASE 1 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 2: Prior-free evaluation at each of the 9 design points
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("PHASE 2: Prior-free evaluation\n")
cat(strrep("═", 70), "\n")

# Assemble design points in density_table order
design_points <- vector("list", 9L)
for (i in seq_along(LEVER_A)) {
  design_points[[i]] <- list(lever="A_globalr", param=sprintf("%.2f",LEVER_A[i]),
                              threshold_r=LEVER_A[i], topk=NULL,
                              edges=edge_lists[[i]])
}
for (j in seq_along(LEVER_B)) {
  design_points[[6L+j]] <- list(lever="B_topk", param=as.character(LEVER_B[j]),
                                 threshold_r=NULL, topk=LEVER_B[j],
                                 edges=edge_lists[[6L+j]])
}

# Check which points already have results (resumability)
already_done <- character(0)
if (file.exists(METRICS_CSV)) {
  prev <- fread(METRICS_CSV)
  if (nrow(prev) > 0L && all(c("lever","param") %in% names(prev))) {
    already_done <- paste0(prev$lever, "_", prev$param)
    stamp(paste("Resuming — skipping", length(already_done), "completed point(s)"))
  }
}

for (pt in design_points) {
  pt_key <- paste0(pt$lever, "_", pt$param)

  if (pt_key %in% already_done) {
    cat("  SKIP", pt_key, "(already in CSV)\n")
    flush.console()
    next
  }

  cat("\n", strrep("-", 60), "\n", sep = "")
  stamp(paste("POINT", pt_key, "| n_pairs =", nrow(pt$edges)))
  t_pt  <- proc.time()
  notes <- character(0)
  timed_out <- FALSE

  # ── 1. Visible genes ──────────────────────────────────────────────────────
  stamp("  computing visible_genes ...")
  vis_res <- .safe(
    stage3_eval_visible_genes(pt$edges, n_total = N_GENES_TOTAL),
    data.frame(n_visible=NA_integer_, n_total=N_GENES_TOTAL, frac_visible=NA_real_),
    "visible_genes"
  )
  cat("    n_visible =", vis_res$n_visible, "\n"); flush.console()

  # ── 2. Per-condition r for retained edges (for fingerprint metrics) ────────
  # Merge from cand_all (all retained edges should be in cand_all)
  edges_with_r <- merge(
    pt$edges[, .(gene_id_A, gene_id_B, abs_r)],
    cand_all[, c("gene_id_A","gene_id_B", COR_COLS), with = FALSE],
    by  = c("gene_id_A","gene_id_B"),
    all.x = TRUE
  )
  n_missing <- sum(is.na(edges_with_r[[COR_COLS[1L]]]))
  if (n_missing > 0L) {
    notes <- c(notes, sprintf("%d edges missing per-cond r", n_missing))
    message("  [INFO] ", n_missing, " edges not in cand_all; fingerprint partial.")
  }

  # ── 3. Effective rank ──────────────────────────────────────────────────────
  stamp("  computing eff_rank ...")
  eff_res <- .safe(
    stage3_eval_effective_rank(edges_with_r, cor_cols = COR_COLS),
    data.frame(eff_rank=NA_real_, n_points=NA_integer_, n_genes=NA_integer_),
    "eff_rank"
  )
  cat("    eff_rank =", round(eff_res$eff_rank, 3), "\n"); flush.console()

  # ── 4. Null gap ────────────────────────────────────────────────────────────
  if (elapsed_secs(t_pt) > BUDGET_SECS) {
    notes <- c(notes, "null_gap skipped: over budget"); timed_out <- TRUE
    cat("  SKIPPING null_gap (over budget)\n"); flush.console()
    null_res <- data.frame(null_gap_ratio=NA_real_, real_frac=NA_real_, perm_frac_mean=NA_real_)
  } else {
    stamp(paste("  computing null_gap (", N_NULL_PERM, "perms) ..."))
    null_res <- .safe(
      stage3_eval_null_gap(edges_with_r, cor_cols=COR_COLS,
                           n_perm=N_NULL_PERM, threshold=0.3, max_genes=2000L),
      data.frame(null_gap_ratio=NA_real_, real_frac=NA_real_, perm_frac_mean=NA_real_),
      "null_gap"
    )
    cat("    null_gap_ratio =", round(null_res$null_gap_ratio, 3), "\n"); flush.console()
  }

  # ── 5. Held-out predictivity ───────────────────────────────────────────────
  if (elapsed_secs(t_pt) > BUDGET_SECS) {
    notes <- c(notes, "heldout_r2 skipped: over budget"); timed_out <- TRUE
    cat("  SKIPPING heldout_predictivity (over budget)\n"); flush.console()
    pred_res <- data.frame(predictivity_mean_r2=NA_real_, predictivity_median_r2=NA_real_)
  } else {
    stamp("  computing heldout_predictivity ...")
    pred_res <- .safe(
      stage3_eval_heldout_predictivity(edges_with_r, cor_cols=COR_COLS, k_partners=10L),
      data.frame(predictivity_mean_r2=NA_real_, predictivity_median_r2=NA_real_),
      "heldout_predictivity"
    )
    cat("    predictivity_mean_r2 =", round(pred_res$predictivity_mean_r2, 4), "\n")
    flush.console()
  }

  # ── 6. Split-half ─────────────────────────────────────────────────────────
  if (elapsed_secs(t_pt) > HARD_TIMEOUT) {
    notes <- c(notes, paste0("TIMED OUT after heldout (", round(elapsed_secs(t_pt)/60,1), " min)"))
    timed_out <- TRUE
    cat("  HARD TIMEOUT — skipping split-half\n"); flush.console()
    sh_res <- data.frame(splithalf_jaccard=NA_real_, splithalf_jaccard_sd=NA_real_, n_splits=0L)
  } else {
    stamp("  computing split-half (3 condition splits) ...")
    # Choose candidate table: Lever A uses abs_r buffer, Lever B uses top-200 subset
    if (!is.null(pt$threshold_r)) {
      cand_sh <- cand_all[abs_r >= (pt$threshold_r - 0.08)]
    } else {
      cand_sh <- cand_topk_sh
    }
    sh_res <- .safe(
      stage3_eval_splithalf(cand_dt=cand_sh, threshold_r=pt$threshold_r,
                            topk=pt$topk, cor_cols=COR_COLS),
      data.frame(splithalf_jaccard=NA_real_, splithalf_jaccard_sd=NA_real_, n_splits=0L),
      "splithalf"
    )
    cat("    splithalf_jaccard =", round(sh_res$splithalf_jaccard, 4),
        "(sd =", round(sh_res$splithalf_jaccard_sd, 4), ")\n")
    flush.console()
  }

  # ── 7. Louvain (descriptive; NOT selection) ────────────────────────────────
  if (elapsed_secs(t_pt) > HARD_TIMEOUT) {
    louvain_res <- data.frame(n_modules=NA_integer_, grey_rate=NA_real_, median_module_size=NA_real_)
  } else {
    stamp("  computing Louvain modules ...")
    louvain_res <- .safe(
      stage3_eval_louvain(pt$edges, seed=98L),
      data.frame(n_modules=NA_integer_, grey_rate=NA_real_, median_module_size=NA_real_),
      "louvain"
    )
    cat("    n_modules =", louvain_res$n_modules,
        "| grey_rate =", round(louvain_res$grey_rate, 3), "\n")
    flush.console()
  }

  eval_secs <- round(elapsed_secs(t_pt), 1)
  cat("  Total elapsed:", round(eval_secs/60, 2), "min\n")

  # ── Write result row ───────────────────────────────────────────────────────
  den_row <- density_table[lever == pt$lever & param == pt$param]
  out_row <- data.frame(
    lever                       = pt$lever,
    param                       = pt$param,
    density                     = den_row$density,
    n_genes                     = vis_res$n_visible,
    splithalf                   = sh_res$splithalf_jaccard,
    eff_rank                    = eff_res$eff_rank,
    null_gap                    = null_res$null_gap_ratio,
    heldout_r2                  = pred_res$predictivity_mean_r2,
    visible_genes               = vis_res$n_visible,
    louvain_n_modules           = louvain_res$n_modules,
    louvain_grey_rate           = louvain_res$grey_rate,
    louvain_median_module_size  = louvain_res$median_module_size,
    eval_seconds                = eval_secs,
    notes                       = paste(notes, collapse = "; "),
    stringsAsFactors            = FALSE
  )
  append_csv_row(out_row, METRICS_CSV)
  stamp(paste("POINT", pt_key, "complete"))
}

cat("\nPHASE 2 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3: Pareto analysis + recommendation
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("PHASE 3: Pareto analysis + recommendation\n")
cat(strrep("═", 70), "\n")

metrics <- fread(METRICS_CSV)

cat("\nMetrics table (key columns):\n")
print(metrics[, .(lever, param, density, splithalf, eff_rank, heldout_r2,
                   visible_genes, louvain_n_modules, louvain_grey_rate)])
flush.console()

# ── 3.1 Pareto front: splithalf (stability) vs heldout_r2 (richness) ────────

is_dominated <- function(dt) {
  sh <- dt$splithalf
  hr <- dt$heldout_r2
  n  <- nrow(dt)
  dom <- logical(n)
  for (i in seq_len(n)) {
    if (is.na(sh[i]) || is.na(hr[i])) { dom[i] <- TRUE; next }
    for (j in seq_len(n)) {
      if (i == j || is.na(sh[j]) || is.na(hr[j])) next
      if (sh[j] >= sh[i] && hr[j] >= hr[i] && (sh[j] > sh[i] || hr[j] > hr[i])) {
        dom[i] <- TRUE; break
      }
    }
  }
  dom
}

metrics[, dominated := is_dominated(.SD)]
pareto_pts <- metrics[dominated == FALSE]

cat("\nPareto-front design points (splithalf vs heldout_r2):\n")
print(pareto_pts[, .(lever, param, density, splithalf, eff_rank, heldout_r2)])
flush.console()

# ── 3.2 Lever A vs Lever B at matched density ────────────────────────────────

cat("\nLever A vs Lever B at matched density:\n")
dens_A <- metrics[lever == "A_globalr"]
dens_B <- metrics[lever == "B_topk"]

lever_compare <- vector("list", nrow(dens_B))
for (bi in seq_len(nrow(dens_B))) {
  d_B   <- dens_B$density[bi]
  k_B   <- dens_B$param[bi]
  diffs <- abs(dens_A$density - d_B)
  best  <- which.min(diffs)
  if (!length(best)) next
  sh_A  <- dens_A$splithalf[best]; sh_B  <- dens_B$splithalf[bi]
  hr_A  <- dens_A$heldout_r2[best]; hr_B  <- dens_B$heldout_r2[bi]
  thr_A <- dens_A$param[best]

  verdict_sh <- if (!is.na(sh_A) && !is.na(sh_B))
                  ifelse(sh_B > sh_A, "B better", ifelse(sh_B < sh_A, "A better", "tie"))
                else "NA"
  verdict_hr <- if (!is.na(hr_A) && !is.na(hr_B))
                  ifelse(hr_B > hr_A, "B better", ifelse(hr_B < hr_A, "A better", "tie"))
                else "NA"

  cat(sprintf("  top-k=%s (density=%.5f) vs |r|>=%.2f (density=%.5f):\n",
              k_B, d_B, as.numeric(thr_A), dens_A$density[best]))
  cat(sprintf("    splithalf:   B=%.4f  A=%.4f  → %s\n", sh_B, sh_A, verdict_sh))
  cat(sprintf("    heldout_r2:  B=%.4f  A=%.4f  → %s\n", hr_B, hr_A, verdict_hr))
  flush.console()

  lever_compare[[bi]] <- data.frame(
    k_B=k_B, d_B=d_B, thr_A=thr_A, d_A=dens_A$density[best],
    sh_A=sh_A, sh_B=sh_B, hr_A=hr_A, hr_B=hr_B,
    verdict_sh=verdict_sh, verdict_hr=verdict_hr,
    stringsAsFactors=FALSE
  )
}

# ── 3.3 Recommendation ───────────────────────────────────────────────────────

if (nrow(pareto_pts) > 0L) {
  # Among Pareto points, prefer highest splithalf; break ties by heldout_r2
  pareto_pts2 <- pareto_pts[!is.na(splithalf) & !is.na(heldout_r2)]
  if (nrow(pareto_pts2) == 0L) pareto_pts2 <- pareto_pts
  # Composite score: splithalf + scaled heldout_r2
  hr_rng <- range(pareto_pts2$heldout_r2, na.rm=TRUE)
  if (diff(hr_rng) > 0)
    pareto_pts2[, composite := splithalf + 0.5 * (heldout_r2 - hr_rng[1]) / diff(hr_rng)]
  else
    pareto_pts2[, composite := splithalf]
  rec <- pareto_pts2[which.max(composite)]
} else {
  rec <- metrics[which.max(splithalf + ifelse(is.na(heldout_r2), 0, heldout_r2 * 0.5))]
}

cat("\n--- RECOMMENDATION (prior-free metrics only) ---\n")
cat("Recommended: lever =", rec$lever, "  param =", rec$param, "\n")
cat("  density        =", round(rec$density, 5), "\n")
cat("  splithalf      =", round(rec$splithalf, 4), " (stability)\n")
cat("  heldout_r2     =", round(rec$heldout_r2, 4), " (richness)\n")
cat("  eff_rank       =", round(rec$eff_rank, 3), "\n")
cat("  visible_genes  =", rec$visible_genes, "\n")
cat("NOTE: This is a recommendation from prior-free metrics only.",
    "The final threshold call is the user's.\n")
flush.console()

# ── 3.4 Post-hoc sanity at recommended point ─────────────────────────────────

cat("\n--- POST-HOC SANITY (at recommended point; NOT a selection input) ---\n")

rec_idx  <- which(sapply(design_points,
                         function(p) paste0(p$lever,"_",p$param)) ==
                    paste0(rec$lever,"_",rec$param))
edges_rec <- if (length(rec_idx)) design_points[[rec_idx[1]]]$edges else
               data.table(gene_id_A=character(), gene_id_B=character(), abs_r=numeric())

net_genes <- unique(c(edges_rec$gene_id_A, edges_rec$gene_id_B))

bon3_partners <- unique(c(
  edges_rec[gene_id_A == BON3_ID, gene_id_B],
  edges_rec[gene_id_B == BON3_ID, gene_id_A]
))
cat("BON3 (AT1G08860):\n")
cat("  in network:", BON3_ID %in% net_genes, "\n")
cat("  n partners:", length(bon3_partners), "\n")

# WRKY genes from symbol map
wrky_ids <- character(0)
sym_paths <- c("results/pathogen_multiome/symbol_map.csv",
               "results/symbol_map.csv")
for (sp in sym_paths) {
  if (file.exists(sp)) {
    sym <- fread(sp)
    # Common column names
    if ("symbol" %in% names(sym) && "gene_id" %in% names(sym)) {
      wrky_ids <- sym[grepl("^WRKY", symbol, ignore.case=FALSE), gene_id]
    } else if (ncol(sym) >= 2) {
      wrky_ids <- sym[[1]][grepl("WRKY", sym[[2]], ignore.case=FALSE)]
    }
    break
  }
}

n_wrky_net <- sum(wrky_ids %in% net_genes)
cat("WRKY family:\n")
cat("  WRKY in universe:", length(wrky_ids), "\n")
cat("  WRKY in network:", n_wrky_net, "\n")
cat("  recovery rate:",
    if (length(wrky_ids) > 0) round(n_wrky_net/length(wrky_ids), 3) else "NA", "\n")
cat("(Post-hoc sanity — NOT used in threshold selection)\n")
flush.console()

# ── 3.5 Write STAGE3_FINDINGS.md ─────────────────────────────────────────────

fmt_dt <- function(dt) {
  cols <- names(dt)
  hdr  <- paste("|", paste(cols, collapse = " | "), "|")
  sep  <- paste("|", paste(rep("---", length(cols)), collapse = " | "), "|")
  rows <- apply(as.data.frame(dt), 1,
                function(r) paste("|", paste(r, collapse = " | "), "|"))
  paste(c(hdr, sep, rows), collapse = "\n")
}

pareto_tbl <- if (nrow(pareto_pts) > 0L)
  fmt_dt(pareto_pts[, .(lever, param, density=round(density,5),
                          splithalf=round(splithalf,4),
                          eff_rank=round(eff_rank,3),
                          heldout_r2=round(heldout_r2,4))])
else "No non-dominated Pareto points (check for NA metrics)."

findings <- c(
  "# Stage 3: Edge-Threshold Selection Findings",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## 1. Density Table (9 design points)",
  "",
  fmt_dt(density_table),
  "",
  "## 2. Prior-Free Metrics Table",
  "",
  fmt_dt(metrics[, .(lever, param, density=round(density,5),
                       splithalf=round(splithalf,4), eff_rank=round(eff_rank,3),
                       null_gap=round(null_gap,2), heldout_r2=round(heldout_r2,4),
                       visible_genes, louvain_n_modules,
                       louvain_grey_rate=round(louvain_grey_rate,3),
                       louvain_median_module_size, eval_seconds, notes)]),
  "",
  "## 3. Pareto Front (stability × richness)",
  "",
  pareto_tbl,
  "",
  "## 4. Recommended Design Point",
  "",
  paste0("**", rec$lever, "  param = ", rec$param, "**"),
  "",
  paste0("| metric | value |"),
  paste0("| --- | --- |"),
  paste0("| density | ", round(rec$density, 5), " |"),
  paste0("| splithalf_jaccard | ", round(rec$splithalf, 4), " |"),
  paste0("| heldout_r2 | ", round(rec$heldout_r2, 4), " |"),
  paste0("| eff_rank | ", round(rec$eff_rank, 3), " |"),
  paste0("| visible_genes | ", rec$visible_genes, " |"),
  "",
  "**Selection rationale:** This design point lies on the stability–richness",
  "Pareto front determined by prior-free metrics only. Split-half Jaccard",
  "(stability) was computed from 3 unique 2+2 condition splits; heldout_r2",
  "from leave-one-condition-out GBA R². No gold-standard gene sets, GO terms,",
  "or motif information were used in the selection.",
  "",
  "> **Final threshold choice is the user's; this is a recommendation",
  "> from prior-free diagnostics only.**",
  "",
  "## 5. Lever A vs Lever B at Matched Density",
  ""
)

for (bi in seq_along(lever_compare)) {
  lc <- lever_compare[[bi]]
  if (is.null(lc)) next
  findings <- c(findings,
    paste0("**top-k=", lc$k_B, " vs |r|>=", lc$thr_A, "** (comparable density)"),
    paste0("- splithalf: B=", round(lc$sh_B,4), "  A=", round(lc$sh_A,4), "  → ", lc$verdict_sh),
    paste0("- heldout_r2: B=", round(lc$hr_B,4), "  A=", round(lc$hr_A,4), "  → ", lc$verdict_hr),
    ""
  )
}

findings <- c(findings,
  "## 6. Post-Hoc Sanity at Recommended Point",
  "",
  "*(NOT used for selection — reported after recommendation is fixed.)*",
  "",
  paste0("**BON3 (AT1G08860):** in network = ", BON3_ID %in% net_genes,
         ";  n partners = ", length(bon3_partners)),
  "",
  paste0("**WRKY family:** ", n_wrky_net, " / ", length(wrky_ids),
         " WRKY genes in network (",
         if (length(wrky_ids) > 0)
           round(100 * n_wrky_net / length(wrky_ids), 1)
         else "NA", "%)"),
  ""
)

writeLines(findings, file.path(OUT_DIR, "STAGE3_FINDINGS.md"))
cat("\nFindings written to", file.path(OUT_DIR, "STAGE3_FINDINGS.md"), "\n")
cat("PHASE 3 complete\n\n")
flush.console()

# ─────────────────────────────────────────────────────────────────────────────
# FINAL REPORT (printed to console)
# ─────────────────────────────────────────────────────────────────────────────

cat(strrep("═", 70), "\n")
cat("FINAL REPORT\n")
cat(strrep("═", 70), "\n\n")

cat("1. DENSITY TABLE (9 rows)\n")
print(density_table)

cat("\n2. METRICS TABLE\n")
print(metrics[, .(lever, param, density, splithalf, eff_rank, heldout_r2,
                   visible_genes, eval_seconds)])

cat("\n3. PARETO FRONT + RECOMMENDATION\n")
if (nrow(pareto_pts) > 0L) {
  print(pareto_pts[, .(lever, param, density, splithalf, heldout_r2)])
} else {
  cat("  (no non-dominated points)\n")
}
cat("Recommended:", rec$lever, " param =", rec$param, "\n")
cat("  stability (splithalf_jaccard):", round(rec$splithalf, 4), "\n")
cat("  richness  (heldout_r2):",        round(rec$heldout_r2, 4), "\n")

cat("\n4. LEVER A vs LEVER B AT MATCHED DENSITY\n")
for (bi in seq_along(lever_compare)) {
  lc <- lever_compare[[bi]]
  if (is.null(lc)) next
  cat(sprintf("  top-k=%s vs |r|>=%.2f: splithalf %s, heldout_r2 %s\n",
              lc$k_B, as.numeric(lc$thr_A), lc$verdict_sh, lc$verdict_hr))
}

cat("\n5. BON3 + WRKY POST-HOC SANITY (NOT selection criteria)\n")
cat("   BON3 in network:", BON3_ID %in% net_genes,
    "| n partners:", length(bon3_partners), "\n")
cat("   WRKY in network:", n_wrky_net, "/", length(wrky_ids), "\n")

cat("\n6. PHASE 4 (metacell sweep): Not attempted — deferred to preserve Phase 3 quality.\n")

cat("\n7. CODE COMMITTED: R/stage3_threshold_eval.R (adapter)\n")
cat("   Results in:", OUT_DIR, "(gitignored)\n")
cat("   Warnings/timeouts: see 'notes' column in", METRICS_CSV, "\n")

cat("\nStage 3 sweep complete.\n")
flush.console()
