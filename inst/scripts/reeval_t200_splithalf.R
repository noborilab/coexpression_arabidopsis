#!/usr/bin/env Rscript
# Re-evaluate metacell_t200 splithalf post-irlba+NaN-fix.
# Patches already in: irlba::prcomp_irlba in .pca_coords (47e7a30)
#                     NaN filter in .matrix_correlation (5dcf50b)

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(irlba)
})

SEURAT_PATH <- file.path(
  path.expand("~"),
  "Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects",
  "SA_PTI_ETI_single_cell",
  "SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed",
  "combined_filtered.rds"
)
SYMBOL_MAP_PATH <- "results/pathogen_multiome/symbol_map.csv"
SUBCLUSTER_COL  <- "sub_clst_rna_20260610"  # dataset-specific: edit for new datasets
CSV_PATH        <- "results/pathogen_multiome/obs_design/stage2_metacell_sweep.csv"

set.seed(98L)
cat(sprintf("[%s] Loading bundle...\n", format(Sys.time())))
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
cat(sprintf("[%s] Bundle: %d genes x %d cells\n", format(Sys.time()),
            nrow(bundle$counts), n_cells))

n_pts <- max(1L, round(n_cells / 200L))
cat(sprintf("[%s] Evaluating metacell_t200 (target_size=200, n_pts=%d)...\n",
            format(Sys.time()), n_pts))
t0 <- proc.time()[["elapsed"]]
set.seed(98L)

result <- evaluate_obs_design(
  bundle        = bundle,
  design_fn     = obs_metacell_knn,
  design_args   = list(target_size = 200L, n_points = n_pts, aggregation = "mean"),
  cor_type      = "spearman",
  norm_method   = "zscore_gene",
  n_splithalf   = 3L,
  heldout_folds = 5L,
  null_perm     = 10L
)

elapsed_s <- round(proc.time()[["elapsed"]] - t0, 1)

cat(sprintf("\n=== t200 RESULT ===\n"))
cat(sprintf("  n_pts             = %d\n",   result$n_points))
cat(sprintf("  eff_rank          = %.3f\n", result$eff_rank))
cat(sprintf("  splithalf_jaccard = %.6f\n", result$splithalf_mat_cor_mean))
cat(sprintf("  heldout_r2        = %.6f\n", result$predictivity_mean_r2))
cat(sprintf("  depth_leakage     = %.6f\n", result$depth_leakage_rho))
cat(sprintf("  visible_genes     = %d\n",   result$n_visible))
cat(sprintf("  eval_seconds      = %.1f\n", elapsed_s))

# Update CSV: replace t200 row
df <- read.csv(CSV_PATH, stringsAsFactors = FALSE)
new_row <- data.frame(
  design            = "metacell_knn",
  target_cells      = 200L,
  n_pts             = result$n_points,
  eff_rank          = result$eff_rank,
  splithalf_jaccard = result$splithalf_mat_cor_mean,
  heldout_r2        = result$predictivity_mean_r2,
  depth_leakage     = result$depth_leakage_rho,
  visible_genes     = result$n_visible,
  eval_seconds      = elapsed_s,
  notes             = "re-eval post-irlba+NaN-fix",
  stringsAsFactors  = FALSE
)
df <- df[df$target_cells != 200L, ]
df <- rbind(new_row, df)
df <- df[order(df$target_cells, decreasing = TRUE), ]
write.csv(df, CSV_PATH, row.names = FALSE)
cat(sprintf("\n[%s] CSV updated: %s\n", format(Sys.time()), CSV_PATH))
cat("REEVAL_T200_DONE\n")
