#!/usr/bin/env Rscript
# FLAG-14 completion: metacell sweep (Arm 1) + cluster granularity re-sweep (Arm 2)
#
# Normalization: zscore_gene + Spearman (FINAL, Stage 1)
# Harness: evaluate_obs_design() from R/coexpr_eval.R (Stage 2 only)
# n_reps=3, seed=98 everywhere
#
# Run:
#   nohup Rscript inst/scripts/flag14_completion_sweep.R \
#     > logs/metacell_sweep.log 2>&1 &
#
# Outputs (not committed — results/ is gitignored):
#   results/pathogen_multiome/obs_design/stage2_metacell_sweep.csv
#   results/pathogen_multiome/obs_design/stage2_cluster_resweep.csv

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)   # picks up irlba patch in .pca_coords
  library(igraph)
  library(irlba)
})

# ---- Constants ----
SEED          <- 98L
NORM_METHOD   <- "zscore_gene"
COR_TYPE      <- "spearman"
N_SPLITHALF   <- 3L
HELDOUT_FOLDS <- 5L
NULL_PERM     <- 10L

TIMEOUT_METACELL <- 45 * 60   # 45 min per design
TIMEOUT_CLUSTER  <- 30 * 60   # 30 min per design

SEURAT_PATH <- file.path(
  path.expand("~"),
  "Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects",
  "SA_PTI_ETI_single_cell",
  "SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed",
  "combined_filtered.rds"
)
SYMBOL_MAP_PATH <- "results/pathogen_multiome/symbol_map.csv"
SUBCLUSTER_COL  <- "sub_clst_rna_20260610"  # dataset-specific: edit for new datasets
OUT_DIR         <- "results/pathogen_multiome/obs_design"

METACELL_CSV <- file.path(OUT_DIR, "stage2_metacell_sweep.csv")
CLUSTER_CSV  <- file.path(OUT_DIR, "stage2_cluster_resweep.csv")

METACELL_SIZES      <- c(200L, 100L, 50L, 25L)
CLUSTER_RESOLUTIONS <- c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0)

# ---- Helpers ----
ts_now <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

append_row_csv <- function(csv_path, row_df) {
  if (file.exists(csv_path)) {
    existing <- read.csv(csv_path, stringsAsFactors = FALSE)
    out_df   <- rbind(existing, row_df)
  } else {
    out_df <- row_df
  }
  write.csv(out_df, csv_path, row.names = FALSE)
}

already_done <- function(csv_path, key_col, key_val) {
  if (!file.exists(csv_path)) return(FALSE)
  existing <- read.csv(csv_path, stringsAsFactors = FALSE)
  isTRUE(key_val %in% existing[[key_col]])
}

# ---- Setup ----
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---- Load bundle ----
cat(sprintf("%s Loading Seurat bundle...\n", ts_now()))
if (!file.exists(SEURAT_PATH)) stop("Seurat object not found: ", SEURAT_PATH)
if (!file.exists(SYMBOL_MAP_PATH)) stop("Symbol map not found: ", SYMBOL_MAP_PATH)

symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)
bundle <- load_seurat(
  seurat_path    = SEURAT_PATH,
  dataset_id     = "pathogen_multiome",
  stratum_var    = "sample2",
  stratum_levels = NULL,
  group_var      = SUBCLUSTER_COL,
  symbol_map     = symbol_map,
  min_cells      = 10L
)
n_cells <- ncol(bundle$counts)
n_genes <- nrow(bundle$counts)
cat(sprintf("%s Bundle: %d genes x %d cells\n", ts_now(), n_genes, n_cells))
if (is.null(bundle$counts_raw)) {
  warning("counts_raw is NULL — depth-downsampling eval not available.")
}


# ============================================================
# ARM 1: Metacell sweep
# obs_metacell_knn with target_size=T, n_points=round(n_cells/T)
# T in c(200, 100, 50, 25)
# ============================================================
cat(sprintf("\n%s === ARM 1: Metacell sweep ===\n", ts_now()))

for (ts in METACELL_SIZES) {

  # Resume: skip if already in CSV
  if (already_done(METACELL_CSV, "target_cells", as.integer(ts))) {
    cat(sprintf("%s METACELL t%d already in CSV — skipping\n", ts_now(), ts))
    next
  }

  n_pts <- max(1L, round(n_cells / as.integer(ts)))
  cat(sprintf("\n%s METACELL t%d  target_size=%d  n_pts~%d\n",
              ts_now(), ts, ts, n_pts))

  t0 <- proc.time()[["elapsed"]]
  set.seed(SEED)

  result <- tryCatch({
    setTimeLimit(elapsed = TIMEOUT_METACELL, transient = TRUE)
    res <- evaluate_obs_design(
      bundle        = bundle,
      design_fn     = obs_metacell_knn,
      design_args   = list(target_size = as.integer(ts),
                           n_points    = n_pts,
                           aggregation = "mean"),
      cor_type      = COR_TYPE,
      norm_method   = NORM_METHOD,
      n_splithalf   = N_SPLITHALF,
      heldout_folds = HELDOUT_FOLDS,
      null_perm     = NULL_PERM
    )
    setTimeLimit(elapsed = Inf)
    res
  }, error = function(e) {
    setTimeLimit(elapsed = Inf)
    msg <- conditionMessage(e)
    if (grepl("time limit|elapsed time", msg, ignore.case = TRUE)) {
      cat(sprintf("%s METACELL t%d TIMED OUT after %.1f min\n",
                  ts_now(), ts, (proc.time()[["elapsed"]] - t0) / 60))
    } else {
      cat(sprintf("%s METACELL t%d FAILED: %s\n", ts_now(), ts, msg))
    }
    NULL
  })

  elapsed_s <- round(proc.time()[["elapsed"]] - t0, 1)

  if (!is.null(result)) {
    row <- data.frame(
      design            = "metacell_knn",
      target_cells      = as.integer(ts),
      n_pts             = result$n_points,
      eff_rank          = result$eff_rank,
      splithalf_jaccard = result$splithalf_mat_cor_mean,
      heldout_r2        = result$predictivity_mean_r2,
      depth_leakage     = result$depth_leakage_rho,
      visible_genes     = result$n_visible,
      eval_seconds      = elapsed_s,
      notes             = "",
      stringsAsFactors  = FALSE
    )
    cat(sprintf("DESIGN metacell_t%d complete — eff_rank=%.1f splithalf=%.3f\n",
                ts, result$eff_rank, result$splithalf_mat_cor_mean))
  } else {
    row <- data.frame(
      design            = "metacell_knn",
      target_cells      = as.integer(ts),
      n_pts             = NA_integer_,
      eff_rank          = NA_real_,
      splithalf_jaccard = NA_real_,
      heldout_r2        = NA_real_,
      depth_leakage     = NA_real_,
      visible_genes     = NA_integer_,
      eval_seconds      = elapsed_s,
      notes             = "TIMED OUT",
      stringsAsFactors  = FALSE
    )
  }

  append_row_csv(METACELL_CSV, row)
  gc(verbose = FALSE)
}


# ============================================================
# ARM 2: Cluster re-sweep
# obs_cluster with resolution R in c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0)
#
# Strategy: pre-compute PCA + kNN graph ONCE for all resolutions, then
# run Louvain per resolution and cache labels in bundle$cell_meta.
# This means split-half reps reuse cached RNA_snn_res.{R} labels and
# skip PCA/kNN entirely — only obs aggregation + correlation remain.
# ============================================================
cat(sprintf("\n%s === ARM 2: Cluster re-sweep ===\n", ts_now()))

resolutions_needed <- CLUSTER_RESOLUTIONS[!sapply(CLUSTER_RESOLUTIONS, function(R) {
  already_done(CLUSTER_CSV, "resolution", R)
})]

if (length(resolutions_needed) > 0L) {
  cat(sprintf("%s Pre-computing PCA + kNN for %d resolutions...\n",
              ts_now(), length(resolutions_needed)))

  # ---- PCA (irlba) on full bundle ----
  n_pcs_cl <- min(30L, n_genes - 1L, n_cells - 1L)
  cat(sprintf("%s   irlba PCA (%d PCs on %d cells)...\n", ts_now(), n_pcs_cl, n_cells))
  pca_cl     <- irlba::prcomp_irlba(t(as.matrix(bundle$counts)), n = n_pcs_cl,
                                     center = TRUE, scale. = FALSE)
  pca_mat_cl <- pca_cl$x   # n_cells x n_pcs
  rm(pca_cl); gc(verbose = FALSE)
  cat(sprintf("%s   PCA done.\n", ts_now()))

  # ---- kNN graph (k=15, chunked) ----
  k_use_cl   <- min(15L, n_cells - 1L)
  norms_cl   <- rowSums(pca_mat_cl^2)
  cs_cl      <- min(500L, n_cells)
  nc_cl      <- ceiling(n_cells / cs_cl)
  fv_cl      <- integer(0L)
  tv_cl      <- integer(0L)

  cat(sprintf("%s   Building kNN graph (%d chunks)...\n", ts_now(), nc_cl))
  for (ci_cl in seq_len(nc_cl)) {
    s_cl <- (ci_cl - 1L) * cs_cl + 1L
    e_cl <- min(ci_cl * cs_cl, n_cells)
    b_cl <- pca_mat_cl[s_cl:e_cl, , drop = FALSE]

    sd_cl <- outer(rowSums(b_cl^2), norms_cl, "+") - 2 * tcrossprod(b_cl, pca_mat_cl)
    sd_cl <- pmax(sd_cl, 0)

    for (li_cl in seq_len(nrow(b_cl))) {
      gi_cl      <- s_cl + li_cl - 1L
      d_cl       <- sd_cl[li_cl, ]
      d_cl[gi_cl] <- Inf
      nn_cl      <- order(d_cl)[seq_len(k_use_cl)]
      fv_cl      <- c(fv_cl, rep(gi_cl, k_use_cl))
      tv_cl      <- c(tv_cl, nn_cl)
    }
  }

  g_cl <- igraph::graph_from_data_frame(
    data.frame(from = fv_cl, to = tv_cl),
    directed = FALSE,
    vertices = data.frame(name = seq_len(n_cells))
  )
  g_cl <- igraph::simplify(g_cl)
  rm(pca_mat_cl, norms_cl, fv_cl, tv_cl, sd_cl, b_cl); gc(verbose = FALSE)
  cat(sprintf("%s   kNN graph done.\n", ts_now()))

  # ---- Louvain per resolution; cache labels in bundle$cell_meta ----
  for (R_pre in resolutions_needed) {
    col_nm <- paste0("RNA_snn_res.", R_pre)
    if (col_nm %in% names(bundle$cell_meta)) next
    cat(sprintf("%s   Louvain res=%.2f...\n", ts_now(), R_pre))
    comm_R <- tryCatch(
      igraph::cluster_louvain(g_cl, resolution = R_pre),
      error = function(e) igraph::cluster_louvain(g_cl)
    )
    bundle$cell_meta[[col_nm]] <- as.character(igraph::membership(comm_R))
    n_cl_R <- length(unique(bundle$cell_meta[[col_nm]]))
    cat(sprintf("%s   res=%.2f -> %d clusters cached in cell_meta\n",
                ts_now(), R_pre, n_cl_R))
    rm(comm_R); gc(verbose = FALSE)
  }
  rm(g_cl); gc(verbose = FALSE)
}

trip_wire_fired <- FALSE

for (R in CLUSTER_RESOLUTIONS) {
  if (trip_wire_fired) break

  # Resume
  if (already_done(CLUSTER_CSV, "resolution", R)) {
    cat(sprintf("%s CLUSTER res=%.2f already in CSV — skipping\n", ts_now(), R))
    next
  }

  cat(sprintf("\n%s CLUSTER res=%.2f\n", ts_now(), R))
  t0 <- proc.time()[["elapsed"]]
  set.seed(SEED)

  result <- tryCatch({
    setTimeLimit(elapsed = TIMEOUT_CLUSTER, transient = TRUE)
    res <- evaluate_obs_design(
      bundle        = bundle,
      design_fn     = obs_cluster,
      design_args   = list(resolution = R, aggregation = "mean"),
      cor_type      = COR_TYPE,
      norm_method   = NORM_METHOD,
      n_splithalf   = N_SPLITHALF,
      heldout_folds = HELDOUT_FOLDS,
      null_perm     = NULL_PERM
    )
    setTimeLimit(elapsed = Inf)
    res
  }, error = function(e) {
    setTimeLimit(elapsed = Inf)
    msg <- conditionMessage(e)
    if (grepl("time limit|elapsed time", msg, ignore.case = TRUE)) {
      cat(sprintf("%s CLUSTER res=%.2f TIMED OUT\n", ts_now(), R))
    } else {
      cat(sprintf("%s CLUSTER res=%.2f FAILED: %s\n", ts_now(), R, msg))
    }
    NULL
  })

  elapsed_s <- round(proc.time()[["elapsed"]] - t0, 1)

  if (!is.null(result)) {
    # Trip wire: n_pts=34 means Bug #1 is still active
    if (R == CLUSTER_RESOLUTIONS[1L] && !is.na(result$n_points) &&
        result$n_points == 34L) {
      cat("BUG #1 STILL ACTIVE — cluster sweep aborted\n")
      trip_wire_fired <- TRUE
    }

    row <- data.frame(
      design            = "obs_cluster",
      resolution        = R,
      n_pts             = result$n_points,
      eff_rank          = result$eff_rank,
      splithalf_jaccard = result$splithalf_mat_cor_mean,
      heldout_r2        = result$predictivity_mean_r2,
      depth_leakage     = result$depth_leakage_rho,
      visible_genes     = result$n_visible,
      eval_seconds      = elapsed_s,
      notes             = if (trip_wire_fired) "BUG#1 DETECTED" else "",
      stringsAsFactors  = FALSE
    )
    cat(sprintf("DESIGN cluster_res%.2f complete — eff_rank=%.1f splithalf=%.3f\n",
                R, result$eff_rank, result$splithalf_mat_cor_mean))
  } else {
    row <- data.frame(
      design            = "obs_cluster",
      resolution        = R,
      n_pts             = NA_integer_,
      eff_rank          = NA_real_,
      splithalf_jaccard = NA_real_,
      heldout_r2        = NA_real_,
      depth_leakage     = NA_real_,
      visible_genes     = NA_integer_,
      eval_seconds      = elapsed_s,
      notes             = "TIMED OUT",
      stringsAsFactors  = FALSE
    )
  }

  append_row_csv(CLUSTER_CSV, row)
  gc(verbose = FALSE)
}

if (trip_wire_fired) {
  cat("\nBUG #1 STILL ACTIVE — cluster sweep aborted after first resolution.\n")
  cat("Proceed to STEP 3 with metacell results only.\n")
}

cat(sprintf("\n%s FLAG-14 completion sweep DONE.\n", ts_now()))
cat(sprintf("  Metacell CSV:  %s\n", METACELL_CSV))
cat(sprintf("  Cluster CSV:   %s\n", CLUSTER_CSV))
