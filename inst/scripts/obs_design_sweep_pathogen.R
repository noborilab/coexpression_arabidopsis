#!/usr/bin/env Rscript
# obs_design_sweep_pathogen.R
#
# Empirical observation-point design exploration on the pathogen multiome data.
# Run under nohup with the OpenBLAS-linked R:
#
#   nohup Rscript inst/scripts/obs_design_sweep_pathogen.R > logs/obs_design_sweep.log 2>&1 &
#
# Stages:
#   1 — Normalization decision: fix one mid-granularity design, test all
#       normalize_obs methods × {spearman, pearson}, pick the default.
#   2 — Granularity sweep: fixed normalization + cor_type, sweep obs_cluster
#       across resolutions and obs_metacell_knn across target_size.
#   Post-hoc sanity (NOT selection): report BON3 + select WRKYs visibility
#       for the 2-3 Pareto-front designs.
#
# Outputs (NOT committed):
#   results/pathogen_multiome/obs_design/normalization_decision.csv
#   results/pathogen_multiome/obs_design/granularity_sweep.csv
#   results/pathogen_multiome/obs_design/OBS_DESIGN_REPORT.md

suppressPackageStartupMessages({
  library(CoexprArabidopsis)   # load_seurat, obs_*, normalize_obs, coexpr_from_obs, eval_*
})

t_start <- proc.time()

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

SEURAT_PATH    <- file.path(
  path.expand("~"),
  "Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects",
  "SA_PTI_ETI_single_cell",
  "SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed",
  "combined_filtered.rds"
)
SYMBOL_MAP_PATH <- "results/pathogen_multiome/symbol_map.csv"
SUBCLUSTER_COL  <- "sub_clst_rna_20260610"  # dataset-specific: edit for new datasets
# sample2 holds 13 time-point×condition values (e.g. "00_Mock", "DC3000_09h").
# Pass NULL so load_seurat auto-detects all levels and keeps all 65k cells.
# Condition context (4-level) is only needed for post-hoc reporting, not design eval.
STRATUM_VAR     <- "sample2"
STRATUM_LEVELS  <- NULL               # keep all 13 sample2 levels
OUT_DIR         <- "results/pathogen_multiome/obs_design"

# Normalization methods to test (Stage 1)
NORM_METHODS <- c("none", "cp10k_log", "log_only", "zscore_gene")
COR_TYPES    <- c("spearman", "pearson")

# Granularity sweep parameters (Stage 2)
CLUSTER_RESOLUTIONS <- c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0)
METACELL_SIZES      <- c(200L, 100L, 50L, 25L)

# Evaluation parameters
N_SPLITHALF_REPS <- 5L
NULL_PERM        <- 20L
HELDOUT_FOLDS    <- 5L
MIN_VAR          <- 1e-6

# Genes of interest for post-hoc sanity (NOT selection criteria)
BON3_ID   <- "AT1G08860"
WRKY_IDS  <- c(WRKY40 = "AT1G80840", WRKY75 = "AT5G13080", WRKY8 = "AT5G46350")

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
cat(sprintf("[%s] obs_design_sweep_pathogen.R starting\n", format(Sys.time())))

if (!file.exists(SEURAT_PATH)) {
  stop("Seurat object not found: ", SEURAT_PATH)
}
if (!file.exists(SYMBOL_MAP_PATH)) {
  stop("Symbol map not found: ", SYMBOL_MAP_PATH)
}

symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)

cat(sprintf("[%s] Loading Seurat object...\n", format(Sys.time())))
bundle <- load_seurat(
  seurat_path    = SEURAT_PATH,
  dataset_id     = "pathogen_multiome",
  stratum_var    = STRATUM_VAR,
  stratum_levels = STRATUM_LEVELS,
  group_var      = SUBCLUSTER_COL,
  symbol_map     = symbol_map,
  min_cells      = 10L
)
cat(sprintf("[%s] Bundle loaded: %d genes x %d cells\n",
            format(Sys.time()), nrow(bundle$counts), ncol(bundle$counts)))
if (is.null(bundle$counts_raw)) {
  warning("counts_raw is NULL — depth-downsampling eval will be skipped.")
}

# ---------------------------------------------------------------------------
# Stage 1 — Normalization decision
# ---------------------------------------------------------------------------
cat(sprintf("[%s] === Stage 1: Normalization decision ===\n", format(Sys.time())))

# Build baseline ObsPointSet using mid-granularity subcluster design
cat(sprintf("[%s]   Building obs_subcluster (col=%s)...\n",
            format(Sys.time()), SUBCLUSTER_COL))
obs_base <- tryCatch(
  obs_subcluster(bundle, group_col = SUBCLUSTER_COL),
  error = function(e) {
    stop("obs_subcluster failed: ", conditionMessage(e))
  }
)
cat(sprintf("[%s]   %d observation points\n", format(Sys.time()), ncol(obs_base$matrix)))

norm_rows <- list()

for (nm in NORM_METHODS) {
  for (ct in COR_TYPES) {
    key <- paste0(nm, "_", ct)
    cat(sprintf("[%s]   Testing norm=%s cor=%s...\n", format(Sys.time()), nm, ct))

    result <- tryCatch({
      # Apply normalization
      obs_norm <- obs_base
      obs_norm$matrix <- normalize_obs(obs_base, method = nm)

      # Eval: depth leakage, split-half, effective rank, visible genes
      dl  <- eval_depth_leakage(obs_norm, threshold = 0.3)
      sh  <- eval_splithalf(bundle,
                             design_fn   = obs_subcluster,
                             design_args = list(group_col = SUBCLUSTER_COL),
                             cor_type    = ct,
                             norm_method = nm,
                             n_reps      = N_SPLITHALF_REPS)
      er  <- eval_effective_rank(obs_norm)
      vg  <- eval_visible_genes(obs_norm, min_var = MIN_VAR)

      data.frame(
        norm_method         = nm,
        cor_type            = ct,
        depth_leakage_rho   = dl$depth_leakage_rho,
        splithalf_mat_cor   = sh$mat_cor_mean,
        splithalf_jaccard   = sh$jaccard_mean,
        eff_rank            = er$eff_rank,
        n_visible           = vg$n_visible,
        frac_visible        = vg$frac_visible,
        stringsAsFactors    = FALSE
      )
    }, error = function(e) {
      warning(sprintf("norm=%s cor=%s failed: %s", nm, ct, conditionMessage(e)))
      data.frame(
        norm_method       = nm,
        cor_type          = ct,
        depth_leakage_rho = NA_real_,
        splithalf_mat_cor = NA_real_,
        splithalf_jaccard = NA_real_,
        eff_rank          = NA_real_,
        n_visible         = NA_integer_,
        frac_visible      = NA_real_,
        stringsAsFactors  = FALSE
      )
    })

    norm_rows[[key]] <- result
    cat(sprintf("[%s]   done (depth_leakage=%.3f, splithalf_cor=%.3f, eff_rank=%.1f)\n",
                format(Sys.time()),
                result$depth_leakage_rho,
                result$splithalf_mat_cor,
                result$eff_rank))
  }
}

norm_df <- do.call(rbind, norm_rows)
rownames(norm_df) <- NULL
norm_csv <- file.path(OUT_DIR, "normalization_decision.csv")
write.csv(norm_df, norm_csv, row.names = FALSE)
cat(sprintf("[%s] Normalization decision table written: %s\n", format(Sys.time()), norm_csv))

# Select default normalization:
# Rule: lowest depth_leakage_rho among methods whose splithalf_mat_cor is within
# 10% of the best (i.e., >= 0.9 * max splithalf_mat_cor).
valid_rows  <- norm_df[!is.na(norm_df$splithalf_mat_cor), ]
best_sh     <- max(valid_rows$splithalf_mat_cor, na.rm = TRUE)
competitive <- valid_rows[valid_rows$splithalf_mat_cor >= 0.9 * best_sh, ]
best_row    <- competitive[which.min(competitive$depth_leakage_rho), ]

DEFAULT_NORM <- best_row$norm_method[[1L]]
DEFAULT_COR  <- best_row$cor_type[[1L]]
cat(sprintf("[%s] CHOSEN NORMALIZATION: norm=%s cor=%s (depth_leakage=%.3f, splithalf=%.3f)\n",
            format(Sys.time()),
            DEFAULT_NORM, DEFAULT_COR,
            best_row$depth_leakage_rho,
            best_row$splithalf_mat_cor))

# ---------------------------------------------------------------------------
# Stage 2 — Granularity sweep
# ---------------------------------------------------------------------------
cat(sprintf("[%s] === Stage 2: Granularity sweep (norm=%s cor=%s) ===\n",
            format(Sys.time()), DEFAULT_NORM, DEFAULT_COR))

sweep_rows <- list()

# Helper: run evaluate_obs_design for one design, catch failures
.run_eval <- function(label, design_fn, design_args) {
  cat(sprintf("[%s]   Evaluating design: %s...\n", format(Sys.time()), label))
  tryCatch({
    res <- evaluate_obs_design(
      bundle          = bundle,
      design_fn       = design_fn,
      design_args     = design_args,
      cor_type        = DEFAULT_COR,
      norm_method     = DEFAULT_NORM,
      n_splithalf     = N_SPLITHALF_REPS,
      splithalf_reps  = N_SPLITHALF_REPS,
      heldout_folds   = HELDOUT_FOLDS,
      null_perm       = NULL_PERM,
      run_downsample_depth = !is.null(bundle$counts_raw),
      run_downsample_cells = TRUE
    )
    res$design_label <- label
    cat(sprintf("[%s]     done (n_points=%d, eff_rank=%.1f, pred_r2=%.3f, splithalf=%.3f)\n",
                format(Sys.time()),
                res$n_points, res$eff_rank,
                res$predictivity_mean_r2, res$splithalf_mat_cor_mean))
    res
  }, error = function(e) {
    warning(sprintf("Design '%s' failed: %s", label, conditionMessage(e)))
    data.frame(design_label = label, design_name = NA_character_,
               n_points = NA_integer_, n_genes = NA_integer_,
               eff_rank = NA_real_, n_visible = NA_integer_,
               frac_visible = NA_real_, predictivity_mean_r2 = NA_real_,
               predictivity_median_r2 = NA_real_, null_gap_ratio = NA_real_,
               real_frac = NA_real_, perm_frac_mean = NA_real_,
               depth_leakage_rho = NA_real_, splithalf_mat_cor_mean = NA_real_,
               splithalf_mat_cor_sd = NA_real_, splithalf_jaccard_mean = NA_real_,
               splithalf_jaccard_sd = NA_real_,
               stringsAsFactors = FALSE)
  })
}

# obs_subcluster baseline
sweep_rows[["subcluster"]] <- .run_eval(
  "subcluster",
  obs_subcluster,
  list(group_col = SUBCLUSTER_COL)
)

# obs_cluster across resolutions
for (res_val in CLUSTER_RESOLUTIONS) {
  label <- sprintf("cluster_res%.2f", res_val)
  sweep_rows[[label]] <- .run_eval(
    label,
    obs_cluster,
    list(resolution = res_val, aggregation = "mean")
  )
}

# obs_metacell_knn across target_size
for (ts in METACELL_SIZES) {
  # n_points: aim for roughly 300 points, capped by ncol(bundle$counts)
  n_pts <- min(300L, ncol(bundle$counts))
  label <- sprintf("metacell_t%d", ts)
  sweep_rows[[label]] <- .run_eval(
    label,
    obs_metacell_knn,
    list(target_size = ts, n_points = n_pts, aggregation = "mean")
  )
}

# Assemble sweep table, align columns
all_cols <- Reduce(union, lapply(sweep_rows, names))
sweep_df <- do.call(rbind, lapply(sweep_rows, function(r) {
  missing <- setdiff(all_cols, names(r))
  r[missing] <- NA
  r[, all_cols, drop = FALSE]
}))
rownames(sweep_df) <- NULL

sweep_csv <- file.path(OUT_DIR, "granularity_sweep.csv")
write.csv(sweep_df, sweep_csv, row.names = FALSE)
cat(sprintf("[%s] Granularity sweep table written: %s\n", format(Sys.time()), sweep_csv))

# ---------------------------------------------------------------------------
# Identify Pareto-front (stability vs richness)
# ---------------------------------------------------------------------------
# Stability: splithalf_mat_cor_mean (higher = better)
# Richness:  predictivity_mean_r2 (higher = better)
# Pareto-front: designs not dominated (no other design has both higher stability
# AND higher richness).

valid_sweep <- sweep_df[!is.na(sweep_df$splithalf_mat_cor_mean) &
                         !is.na(sweep_df$predictivity_mean_r2), ]

is_pareto <- function(sh, pred) {
  n <- length(sh)
  pareto <- logical(n)
  for (i in seq_len(n)) {
    dominated <- any(sh[-i] >= sh[i] & pred[-i] >= pred[i] &
                       (sh[-i] > sh[i] | pred[-i] > pred[i]))
    pareto[i] <- !dominated
  }
  pareto
}

if (nrow(valid_sweep) > 0L) {
  valid_sweep$is_pareto <- is_pareto(
    valid_sweep$splithalf_mat_cor_mean,
    valid_sweep$predictivity_mean_r2
  )
  pareto_designs <- valid_sweep[valid_sweep$is_pareto, ]
  cat(sprintf("[%s] Pareto-front designs (%d):\n", format(Sys.time()), nrow(pareto_designs)))
  for (i in seq_len(nrow(pareto_designs))) {
    cat(sprintf("  %s: splithalf=%.3f, pred_r2=%.3f, eff_rank=%.1f\n",
                pareto_designs$design_label[i],
                pareto_designs$splithalf_mat_cor_mean[i],
                pareto_designs$predictivity_mean_r2[i],
                pareto_designs$eff_rank[i]))
  }
} else {
  pareto_designs <- valid_sweep
  valid_sweep$is_pareto <- FALSE
}

# Pick top-3 Pareto-front designs for post-hoc sanity
top_pareto <- head(pareto_designs[order(-pareto_designs$predictivity_mean_r2), ], 3)

# ---------------------------------------------------------------------------
# Post-hoc sanity readout (NOT selection criteria)
# ---------------------------------------------------------------------------
cat(sprintf("[%s] === Post-hoc sanity (NOT selection) ===\n", format(Sys.time())))

goi_ids   <- c(BON3 = BON3_ID, WRKY_IDS)
goi_report <- list()

for (i in seq_len(nrow(top_pareto))) {
  label <- top_pareto$design_label[i]
  cat(sprintf("[%s]   Post-hoc for design: %s\n", format(Sys.time()), label))

  # Reconstruct the ObsPointSet for this design
  obs_posthoc <- tryCatch({
    parts  <- strsplit(label, "_")[[1L]]
    prefix <- paste(parts[1:2], collapse="_")

    obs_raw <- if (startsWith(label, "subcluster")) {
      obs_subcluster(bundle, group_col = SUBCLUSTER_COL)
    } else if (startsWith(label, "cluster_res")) {
      res_val <- as.numeric(sub("cluster_res", "", label))
      obs_cluster(bundle, resolution = res_val, aggregation = "mean")
    } else if (startsWith(label, "metacell_t")) {
      ts <- as.integer(sub("metacell_t", "", label))
      obs_metacell_knn(bundle, target_size = ts, n_points = min(300L, ncol(bundle$counts)))
    } else {
      NULL
    }

    if (!is.null(obs_raw)) {
      obs_raw$matrix <- normalize_obs(obs_raw, method = DEFAULT_NORM)
    }
    obs_raw
  }, error = function(e) {
    warning("Post-hoc obs build failed for ", label, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(obs_posthoc)) next

  # Compute correlations
  coexpr <- tryCatch(
    coexpr_from_obs(obs_posthoc, cor_type = DEFAULT_COR, storage_cutoff = 0.1),
    error = function(e) NULL
  )
  if (is.null(coexpr)) next

  et <- coexpr$edge_table

  for (goi_name in names(goi_ids)) {
    gid <- goi_ids[[goi_name]]
    visible <- gid %in% obs_posthoc$gene_ids
    var_val <- if (visible) var(obs_posthoc$matrix[gid, ]) else NA_real_
    n_partners <- if (visible && nrow(et) > 0L) {
      sum(et$gene_id_A == gid | et$gene_id_B == gid)
    } else {
      NA_integer_
    }
    # Top-3 partners
    if (visible && nrow(et) > 0L) {
      partners_et <- et[(et$gene_id_A == gid | et$gene_id_B == gid), ]
      partners_et <- partners_et[order(-abs(partners_et$weight)), ]
      top3 <- head(partners_et, 3L)
      top3_str <- paste(
        ifelse(top3$gene_id_A == gid, top3$gene_id_B, top3$gene_id_A),
        sprintf("(%.3f)", top3$weight),
        sep = " ",
        collapse = "; "
      )
    } else {
      top3_str <- NA_character_
    }

    goi_report[[paste(label, goi_name, sep = "|")]] <- data.frame(
      design      = label,
      gene        = goi_name,
      gene_id     = gid,
      visible     = visible,
      variance    = var_val,
      n_partners  = n_partners,
      top3_partners = top3_str,
      stringsAsFactors = FALSE
    )
  }
}

goi_df <- if (length(goi_report) > 0L) {
  do.call(rbind, goi_report)
} else {
  data.frame(design=character(), gene=character(), gene_id=character(),
             visible=logical(), variance=numeric(), n_partners=integer(),
             top3_partners=character(), stringsAsFactors=FALSE)
}
rownames(goi_df) <- NULL

# ---------------------------------------------------------------------------
# Write OBS_DESIGN_REPORT.md
# ---------------------------------------------------------------------------
t_elapsed <- (proc.time() - t_start)[["elapsed"]]
total_min  <- round(t_elapsed / 60, 1)

report_lines <- c(
  "# Observation-Point Design Sweep Report — Pathogen Multiome",
  "",
  sprintf("Generated: %s | Total wall time: %.1f min", format(Sys.time()), total_min),
  sprintf("Seurat object: %s", SEURAT_PATH),
  sprintf("Subcluster column: %s", SUBCLUSTER_COL),
  "",
  "---",
  "",
  "## Stage 1: Normalization Decision",
  "",
  sprintf("Selection rule: lowest depth_leakage_rho among methods with splithalf_mat_cor >= %.0f%% of best.",
          90),
  sprintf("Chosen: norm=%s, cor=%s", DEFAULT_NORM, DEFAULT_COR),
  "",
  "Full table (all normalization × cor_type combinations):",
  "",
  "| norm_method | cor_type | depth_leakage_rho | splithalf_mat_cor | eff_rank | n_visible |",
  "|---|---|---|---|---|---|",
  apply(norm_df, 1, function(r) {
    sprintf("| %s | %s | %.3f | %.3f | %.1f | %s |",
            r["norm_method"], r["cor_type"],
            as.numeric(r["depth_leakage_rho"]),
            as.numeric(r["splithalf_mat_cor"]),
            as.numeric(r["eff_rank"]),
            r["n_visible"])
  }),
  "",
  "---",
  "",
  "## Stage 2: Granularity Sweep",
  "",
  sprintf("Normalization: %s | Correlation: %s", DEFAULT_NORM, DEFAULT_COR),
  "",
  "| design | n_points | eff_rank | n_visible | pred_r2 | splithalf_cor | depth_leak | is_pareto |",
  "|---|---|---|---|---|---|---|---|"
)

if (nrow(valid_sweep) > 0L) {
  sweep_for_report <- merge(sweep_df, valid_sweep[, c("design_label", "is_pareto")],
                             by = "design_label", all.x = TRUE)
  sweep_for_report$is_pareto[is.na(sweep_for_report$is_pareto)] <- FALSE
  for (i in seq_len(nrow(sweep_for_report))) {
    r <- sweep_for_report[i, ]
    report_lines <- c(report_lines,
      sprintf("| %s | %s | %.1f | %s | %.3f | %.3f | %.3f | %s |",
              r$design_label,
              ifelse(is.na(r$n_points), "NA", as.character(r$n_points)),
              ifelse(is.na(r$eff_rank), NA_real_, r$eff_rank),
              ifelse(is.na(r$n_visible), "NA", as.character(r$n_visible)),
              ifelse(is.na(r$predictivity_mean_r2), NA_real_, r$predictivity_mean_r2),
              ifelse(is.na(r$splithalf_mat_cor_mean), NA_real_, r$splithalf_mat_cor_mean),
              ifelse(is.na(r$depth_leakage_rho), NA_real_, r$depth_leakage_rho),
              ifelse(isTRUE(r$is_pareto), "**YES**", "no")))
  }
}

report_lines <- c(report_lines,
  "",
  "### Stability-Richness Pareto Front",
  "",
  "Stability = splithalf_mat_cor_mean (reproducibility); Richness = predictivity_mean_r2.",
  "Pareto-front designs are not dominated on both axes simultaneously.",
  ""
)

if (nrow(top_pareto) > 0L) {
  for (i in seq_len(nrow(top_pareto))) {
    r <- top_pareto[i, ]
    report_lines <- c(report_lines,
      sprintf("- **%s**: splithalf=%.3f, pred_r2=%.3f, eff_rank=%.1f, n_points=%s",
              r$design_label, r$splithalf_mat_cor_mean, r$predictivity_mean_r2,
              r$eff_rank, r$n_points))
  }
} else {
  report_lines <- c(report_lines, "(No valid designs to compare.)")
}

report_lines <- c(report_lines,
  "",
  "---",
  "",
  "## Post-Hoc Sanity Readout",
  "",
  "> **WARNING**: This section is a sanity check ONLY. BON3 and WRKY visibility",
  "> are NOT selection criteria. Design was chosen on prior-free metrics above.",
  "",
  "Genes checked: BON3 (AT1G08860), WRKY40 (AT1G80840), WRKY75 (AT5G13080), WRKY8 (AT5G46350).",
  "",
  "| design | gene | gene_id | visible | variance | n_partners | top3_partners |",
  "|---|---|---|---|---|---|---|"
)

if (nrow(goi_df) > 0L) {
  for (i in seq_len(nrow(goi_df))) {
    r <- goi_df[i, ]
    report_lines <- c(report_lines,
      sprintf("| %s | %s | %s | %s | %.4f | %s | %s |",
              r$design, r$gene, r$gene_id,
              ifelse(isTRUE(r$visible), "YES", "no"),
              ifelse(is.na(r$variance), 0, r$variance),
              ifelse(is.na(r$n_partners), "NA", as.character(r$n_partners)),
              ifelse(is.na(r$top3_partners), "NA", r$top3_partners)))
  }
}

report_lines <- c(report_lines,
  "",
  "---",
  "",
  "## Notes",
  "",
  "- All design selection decisions are based solely on prior-free metrics",
  "  (stability, richness, null gap, depth leakage).",
  "- The post-hoc BON3/WRKY readout above is a sanity check to confirm that",
  "  the chosen design(s) resolve genes known to have biologically meaningful",
  "  co-expression. It is NOT part of the selection criterion.",
  "- Normalization is an open empirical question (FLAG-14); this sweep settles",
  "  it for this dataset and should be re-run for new datasets.",
  sprintf("- Results produced with CoexprArabidopsis (wall time: %.1f min).", total_min)
)

report_path <- file.path(OUT_DIR, "OBS_DESIGN_REPORT.md")
writeLines(report_lines, report_path)
cat(sprintf("[%s] Report written: %s\n", format(Sys.time()), report_path))

cat(sprintf("\n[%s] obs_design_sweep_pathogen.R COMPLETE (%.1f min)\n",
            format(Sys.time()), total_min))
cat(sprintf("  normalization_decision.csv : %s\n", norm_csv))
cat(sprintf("  granularity_sweep.csv      : %s\n", sweep_csv))
cat(sprintf("  OBS_DESIGN_REPORT.md       : %s\n", report_path))
