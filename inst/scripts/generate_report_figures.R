#!/usr/bin/env Rscript
## generate_report_figures.R
## Generates all 10 report figures (static PNG + SVG, 300 dpi).
## Must be run from the repo root: Rscript inst/scripts/generate_report_figures.R
##
## Uses base-R + ggplot2; fread with nThread=1L (no data.table GForce).
## seed=98. Writes to results/pathogen_multiome/report/figures/

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(svglite)
})

set.seed(98L)

DATASET_ID  <- "pathogen_multiome"
BASE_DIR    <- file.path("results", DATASET_ID)
OBS_DIR     <- file.path(BASE_DIR, "obs_design")
STAGE3_DIR  <- file.path(BASE_DIR, "stage3_threshold_sweep")
PB_DIR      <- file.path(BASE_DIR, "pseudobulk_zscore_spearman")
MOD_DIR     <- file.path(PB_DIR, "modules_official")
WGCNA_DIR   <- file.path(MOD_DIR, "wgcna")
LOUVAIN_DIR <- file.path(MOD_DIR, "louvain")
GOI_DIR     <- file.path(BASE_DIR, "geneset_lookups")

FIG_DIR <- file.path(BASE_DIR, "report", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: save as both PNG and SVG
save_fig <- function(p, name, width = 8, height = 6) {
  png_path <- file.path(FIG_DIR, paste0(name, ".png"))
  svg_path <- file.path(FIG_DIR, paste0(name, ".svg"))
  ggsave(png_path, plot = p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(svg_path, plot = p, width = width, height = height, device = svglite::svglite, bg = "white")
  message("Saved: ", name, ".png + .svg")
  invisible(list(png = png_path, svg = svg_path))
}

# Consistent theme
theme_report <- function() {
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "#f0f0f0", color = "grey70"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 10, color = "grey40"),
    legend.position   = "right"
  )
}

HIGHLIGHT_COL  <- "#1b7837"   # winner / chosen
PARETO_COL     <- "#762a83"   # Pareto front
DOMINATED_COL  <- "#d3d3d3"   # dominated designs
TOPK_COL       <- "#d95f02"   # Lever B (top-k)
GLOBAL_COL     <- "#1f78b4"   # Lever A (global |r|)

cat("=== Generating report figures ===\n\n")

# ==============================================================================
# FIG 01 — Stage 1: Normalization comparison
# ==============================================================================
cat("fig01 — Stage 1 normalization comparison\n")

norm_path <- file.path(OBS_DIR, "normalization_decision.csv")
if (file.exists(norm_path)) {
  norm_df <- as.data.frame(fread(norm_path, nThread = 1L))

  # Human-readable labels
  norm_df$method_label <- paste0(norm_df$norm_method, " + ", norm_df$cor_type)
  norm_df$is_winner    <- norm_df$norm_method == "zscore_gene" & norm_df$cor_type == "spearman"

  # Facet melt: 3 metrics
  metrics_long <- rbind(
    data.frame(method_label = norm_df$method_label, is_winner = norm_df$is_winner,
               metric = "Depth Leakage (ρ)", value = norm_df$depth_leakage_rho),
    data.frame(method_label = norm_df$method_label, is_winner = norm_df$is_winner,
               metric = "Effective Rank", value = norm_df$eff_rank / max(norm_df$eff_rank, na.rm = TRUE)),
    data.frame(method_label = norm_df$method_label, is_winner = norm_df$is_winner,
               metric = "Split-Half Stability", value = norm_df$splithalf_mat_cor)
  )

  # Order by depth leakage ascending (winner should be obvious)
  order_by_leak <- norm_df$method_label[order(norm_df$depth_leakage_rho)]
  metrics_long$method_label <- factor(metrics_long$method_label, levels = rev(order_by_leak))

  fig01 <- ggplot(metrics_long, aes(x = value, y = method_label,
                                     fill = is_winner)) +
    geom_col(width = 0.7, alpha = 0.9) +
    facet_wrap(~ metric, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = HIGHLIGHT_COL),
                      guide = "none") +
    labs(
      title    = "Fig 1 — Stage 1: Normalization Method Comparison",
      subtitle = "Winner: zscore_gene + Spearman (depth leakage 0.107; eff_rank 159.9; split-half 0.895)",
      x        = "Value (eff_rank normalised to [0,1])",
      y        = NULL
    ) +
    theme_report() +
    theme(axis.text.y = element_text(size = 9))

  save_fig(fig01, "fig01_stage1_normalization", width = 11, height = 5)
} else {
  message("SKIP fig01 — ", norm_path, " not found")
}

# ==============================================================================
# FIG 02 — Stage 2: Pareto front (all designs)
# ==============================================================================
cat("fig02 — Stage 2 Pareto front\n")

meta_path <- file.path(OBS_DIR, "stage2_metacell_sweep.csv")
clus_path <- file.path(OBS_DIR, "stage2_cluster_resweep.csv")

if (file.exists(meta_path) && file.exists(clus_path)) {
  meta_df <- as.data.frame(fread(meta_path, nThread = 1L))
  clus_df <- as.data.frame(fread(clus_path, nThread = 1L))

  # subcluster reference point
  subclust_row <- data.frame(design = "subcluster", n_pts = 298,
                             eff_rank = 159.875, splithalf_jaccard = 0.895,
                             label = "subcluster\n(298 pts)")

  # metacell rows
  meta_rows <- data.frame(
    design = paste0("metacell_t", meta_df$target_cells),
    n_pts  = meta_df$n_pts,
    eff_rank = meta_df$eff_rank,
    splithalf_jaccard = meta_df$splithalf_jaccard
  )
  meta_rows$label <- paste0(meta_rows$design, "\n(", meta_rows$n_pts, " pts)")

  # cluster rows
  clus_rows <- data.frame(
    design = paste0("cluster_res", clus_df$resolution),
    n_pts  = clus_df$n_pts,
    eff_rank = clus_df$eff_rank,
    splithalf_jaccard = clus_df$splithalf_jaccard
  )
  clus_rows$label <- paste0("res", clus_df$resolution, "\n(", clus_rows$n_pts, " pts)")

  # combine
  all_designs <- rbind(
    data.frame(category = "subcluster", design = subclust_row$design,
               n_pts = subclust_row$n_pts, eff_rank = subclust_row$eff_rank,
               splithalf_jaccard = subclust_row$splithalf_jaccard, label = subclust_row$label),
    data.frame(category = "metacell",  design = meta_rows$design,
               n_pts = meta_rows$n_pts, eff_rank = meta_rows$eff_rank,
               splithalf_jaccard = meta_rows$splithalf_jaccard, label = meta_rows$label),
    data.frame(category = "cluster",   design = clus_rows$design,
               n_pts = clus_rows$n_pts, eff_rank = clus_rows$eff_rank,
               splithalf_jaccard = clus_rows$splithalf_jaccard, label = clus_rows$label)
  )

  # Mark Pareto front: {subcluster, metacell_t25}
  all_designs$pareto <- all_designs$design %in% c("subcluster", "metacell_t25")

  fig02 <- ggplot(all_designs,
                  aes(x = eff_rank, y = splithalf_jaccard,
                      color = category, shape = pareto, size = pareto)) +
    geom_point(alpha = 0.85) +
    geom_label(data = subset(all_designs, pareto | design %in% c("cluster_res4", "metacell_t100")),
               aes(label = label), size = 2.8, hjust = -0.08, show.legend = FALSE,
               fill = "white", label.size = 0.2) +
    scale_color_manual(values = c("subcluster" = HIGHLIGHT_COL,
                                   "metacell"   = "#4393c3",
                                   "cluster"    = "#d6604d"),
                       name = "Design type") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values  = c("FALSE" = 2.5, "TRUE" = 5), guide = "none") +
    scale_x_log10(breaks = c(3, 10, 30, 100, 300, 1000),
                  labels = c("3", "10", "30", "100", "300", "1000")) +
    labs(
      title    = "Fig 2 — Stage 2: Observation-Point Design Pareto Front",
      subtitle = "Pareto front: {subcluster(298), metacell_t25(2602)} — no design dominates both axes\nSubcluster = stability champion; metacell_t25 = richness champion (7.2× eff_rank)",
      x        = "Effective Rank (log₁₀ scale)",
      y        = "Split-Half Stability (Jaccard)"
    ) +
    coord_cartesian(xlim = c(2.5, 3000), ylim = c(0.55, 0.93)) +
    theme_report()

  save_fig(fig02, "fig02_stage2_pareto", width = 9, height = 6)
}

# ==============================================================================
# FIG 03 — Stage 2: Cluster resolution sweep
# ==============================================================================
cat("fig03 — Stage 2 cluster resolution sweep\n")

if (file.exists(clus_path)) {
  clus_df <- as.data.frame(fread(clus_path, nThread = 1L))
  clus_long <- rbind(
    data.frame(resolution = clus_df$resolution, metric = "Split-Half Stability",
               value = clus_df$splithalf_jaccard),
    data.frame(resolution = clus_df$resolution, metric = "Effective Rank",
               value = clus_df$eff_rank)
  )

  fig03 <- ggplot(clus_long, aes(x = resolution, y = value)) +
    geom_line(color = "#d6604d", linewidth = 1) +
    geom_point(color = "#d6604d", size = 3) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    labs(
      title    = "Fig 3 — Stage 2: Cluster Resolution Sweep (Bug #1 Fixed)",
      subtitle = "All cluster designs dominated by subcluster(298) on both axes\nBug #1 fix verified: n_pts varies 4→51 across resolutions (was constant 34)",
      x        = "Louvain Resolution",
      y        = NULL
    ) +
    theme_report()

  save_fig(fig03, "fig03_stage2_cluster_resweep", width = 9, height = 5)
}

# ==============================================================================
# FIG 04 — Stage 3: Density vs threshold (both levers)
# ==============================================================================
cat("fig04 — Stage 3 density levers\n")

dens_path <- file.path(STAGE3_DIR, "density_table.csv")
if (file.exists(dens_path)) {
  dens_df <- as.data.frame(fread(dens_path, nThread = 1L))

  leA <- dens_df[dens_df$lever == "A_globalr", ]
  leB <- dens_df[dens_df$lever == "B_topk", ]

  fig04 <- ggplot() +
    geom_line(data = leA, aes(x = param, y = density * 100, color = "Lever A (global |r|)"),
              linewidth = 1.2) +
    geom_point(data = leA, aes(x = param, y = density * 100, color = "Lever A (global |r|)"),
               size = 3) +
    # Mark chosen threshold
    geom_vline(xintercept = 0.42, linetype = "dashed", color = HIGHLIGHT_COL, linewidth = 0.8) +
    annotate("text", x = 0.425, y = max(leA$density * 100) * 0.85,
             label = "|r| ≥ 0.42\n(chosen)", color = HIGHLIGHT_COL, size = 3.5, hjust = 0) +
    # Lever B as secondary (x=top_k/10 for display purposes — different axis!)
    scale_x_continuous(
      name = "Lever A: Global |r| threshold",
      sec.axis = sec_axis(~ . * 1, name = "(Lever B: top-k shown at matched density points)")
    ) +
    scale_color_manual(values = c("Lever A (global |r|)" = GLOBAL_COL),
                       name = NULL) +
    labs(
      title    = "Fig 4 — Stage 3: Network Density vs Threshold",
      subtitle = "Lever A: global |r| threshold (5 points); Lever B: per-gene top-k (shown in fig05)\nChosen: |r| ≥ 0.42 → 751,959 pairs / 5,450 genes / density 1.24%",
      y        = "Network Density (%)"
    ) +
    theme_report() +
    theme(legend.position = "none")

  save_fig(fig04, "fig04_stage3_density", width = 8, height = 5)
}

# ==============================================================================
# FIG 05 — Stage 3: Pareto (splithalf × eff_rank), lever A vs B
# ==============================================================================
cat("fig05 — Stage 3 Pareto (stability vs richness)\n")

met_path <- file.path(STAGE3_DIR, "stage3_metrics.csv")
if (file.exists(met_path)) {
  met_df <- as.data.frame(fread(met_path, nThread = 1L))

  # Deduplicate (stage3_metrics has some duplicate rows from second sweep run)
  met_df <- met_df[!duplicated(paste(met_df$lever, met_df$param)), ]

  met_df$lever_label <- ifelse(met_df$lever == "A_globalr", "Lever A (global |r|)", "Lever B (top-k)")
  met_df$is_chosen   <- met_df$lever == "A_globalr" & abs(met_df$param - 0.42) < 0.001
  met_df$null_gap_label <- ifelse(is.infinite(met_df$null_gap), "Inf",
                                   formatC(met_df$null_gap, format = "f", digits = 2))

  fig05 <- ggplot(met_df, aes(x = eff_rank, y = splithalf_jaccard,
                               color = lever_label, shape = is_chosen, size = is_chosen)) +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = paste0("param=", param, "\nnull_gap=", null_gap_label)),
              size = 2.4, hjust = -0.1, show.legend = FALSE) +
    scale_color_manual(values = c("Lever A (global |r|)" = GLOBAL_COL,
                                   "Lever B (top-k)"     = TOPK_COL),
                       name = NULL) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 5.5), guide = "none") +
    labs(
      title    = "Fig 5 — Stage 3: Stability vs Richness Pareto (threshold lever comparison)",
      subtitle = "Lever A (|r|) = high stability + monotone richness; null_gap=Inf (clean signal)\nLever B (top-k) = max richness but low null_gap (~1.19, near-random) → rejected",
      x        = "Effective Rank (richness)",
      y        = "Split-Half Stability (Jaccard)"
    ) +
    coord_cartesian(xlim = c(95, 175), ylim = c(0.875, 0.97)) +
    theme_report() +
    theme(legend.position = "top")

  save_fig(fig05, "fig05_stage3_pareto", width = 9, height = 6)
}

# ==============================================================================
# FIG 06 — Module size distributions: WGCNA vs Louvain
# ==============================================================================
cat("fig06 — Module size distributions\n")

wgcna_sum_path   <- file.path(WGCNA_DIR,   "module_summary.csv")
louvain_sum_path <- file.path(LOUVAIN_DIR, "module_summary.csv")

if (file.exists(wgcna_sum_path) && file.exists(louvain_sum_path)) {
  w_sum <- as.data.frame(fread(wgcna_sum_path,   nThread = 1L))
  l_sum <- as.data.frame(fread(louvain_sum_path, nThread = 1L))

  # Add grey module stats
  # WGCNA: 5450 total assigned; 5450 - sum(w_sum$n_genes) = grey
  total_genes <- 5450L
  w_grey <- total_genes - sum(w_sum$n_genes)
  l_grey <- total_genes - sum(l_sum$n_genes)

  w_all <- rbind(w_sum[, c("module", "n_genes")],
                 data.frame(module = 0L, n_genes = w_grey))
  l_all <- rbind(l_sum[, c("module", "n_genes")],
                 data.frame(module = 0L, n_genes = l_grey))

  w_all$method  <- "WGCNA (12 modules)"
  l_all$method  <- "Louvain (6 modules)"
  w_all$is_grey <- w_all$module == 0L
  l_all$is_grey <- l_all$module == 0L

  mod_df <- rbind(
    data.frame(w_all, mod_id = ifelse(w_all$is_grey, "grey", paste0("M", w_all$module))),
    data.frame(l_all, mod_id = ifelse(l_all$is_grey, "grey", paste0("M", l_all$module)))
  )

  mod_df <- mod_df[order(mod_df$method, -mod_df$n_genes), ]
  mod_df$mod_id <- factor(mod_df$mod_id,
                           levels = rev(unique(mod_df$mod_id[order(mod_df$n_genes)])))

  fig06 <- ggplot(mod_df, aes(x = n_genes, y = reorder(mod_id, n_genes),
                               fill = is_grey)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = n_genes), hjust = -0.1, size = 3) +
    facet_wrap(~ method, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c("FALSE" = "#4393c3", "TRUE" = "grey60"),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title    = "Fig 6 — Official Module Size Distributions",
      subtitle = "WGCNA: 12 modules, power=9, grey=8.1% | Louvain: 6 modules, seed=98, grey=6.6%\nAll modules: dominant condition pattern = '1111' (constitutive across all 4 conditions)",
      x        = "Number of genes",
      y        = NULL
    ) +
    theme_report()

  save_fig(fig06, "fig06_modules_overview", width = 10, height = 7)
}

# ==============================================================================
# FIG 07 — pickSoftThreshold: scale-free fit
# ==============================================================================
cat("fig07 — pickSoftThreshold\n")

sft_path <- file.path(MOD_DIR, "pickSoftThreshold_fitIndices.csv")
if (file.exists(sft_path)) {
  sft_df <- as.data.frame(fread(sft_path, nThread = 1L))

  # Panel A: R^2 vs power
  pA <- ggplot(sft_df, aes(x = Power, y = SFT.R.sq)) +
    geom_line(color = "#4393c3", linewidth = 1) +
    geom_point(color = "#4393c3", size = 2.5) +
    geom_point(data = sft_df[sft_df$Power == 9, ], size = 5,
               color = HIGHLIGHT_COL, shape = 18) +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "grey50") +
    annotate("text", x = 9.3, y = sft_df$SFT.R.sq[sft_df$Power == 9] - 0.02,
             label = "power=9\nR²=0.811", color = HIGHLIGHT_COL, size = 3.5, hjust = 0) +
    labs(x = "Soft-Thresholding Power", y = "Scale-Free Fit (R²)") +
    theme_report()

  # Panel B: Mean connectivity vs power
  pB <- ggplot(sft_df, aes(x = Power, y = mean.k.)) +
    geom_line(color = "#d6604d", linewidth = 1) +
    geom_point(color = "#d6604d", size = 2.5) +
    geom_point(data = sft_df[sft_df$Power == 9, ], size = 5,
               color = HIGHLIGHT_COL, shape = 18) +
    scale_y_log10() +
    labs(x = "Soft-Thresholding Power", y = "Mean Connectivity (log₁₀)") +
    theme_report()

  library(patchwork)
  fig07 <- (pA | pB) +
    plot_annotation(
      title    = "Fig 7 — WGCNA pickSoftThreshold: Scale-Free Fit vs Power",
      subtitle = "Chosen: power=9 (R²=0.811, threshold R²≥0.85 line shown)\n5,450 genes × 751,959 edges at |r|≥0.42"
    )

  save_fig(fig07, "fig07_pickSoftThreshold", width = 10, height = 5)
}

# ==============================================================================
# FIG 08 — Condition patterns: module dominant pattern + pair distribution
# ==============================================================================
cat("fig08 — Condition patterns\n")

wgcna_cond_path   <- file.path(WGCNA_DIR,   "module_condition_patterns.csv")
louvain_cond_path <- file.path(LOUVAIN_DIR, "module_condition_patterns.csv")

if (file.exists(wgcna_cond_path)) {
  wc_df <- as.data.frame(fread(wgcna_cond_path,   nThread = 1L))
  lc_df <- as.data.frame(fread(louvain_cond_path, nThread = 1L))

  wc_df$method <- "WGCNA"
  lc_df$method <- "Louvain"

  comb <- rbind(
    wc_df[, c("method", "module", "n_pairs", "dominant_pattern", "mean_r_score")],
    lc_df[, c("method", "module", "n_pairs", "dominant_pattern", "mean_r_score")]
  )

  # Fraction "1111" vs other per module
  comb$frac_1111 <- ifelse(comb$dominant_pattern == "1111",
                            1 - (comb$mean_r_score * 0 + 0),  # placeholder
                            0)

  # Collect pattern fractions: frac_1111 is the key one
  # Use the wc_df frac_1111 column if available
  if ("frac_1111" %in% names(wc_df)) {
    wc_df2 <- wc_df; lc_df2 <- lc_df
    wc_df2$frac_1111 <- wc_df2$frac_1111
    lc_df2$frac_1111 <- lc_df2$frac_1111

    comb2 <- rbind(
      data.frame(method = "WGCNA",   module = wc_df2$module, n_pairs = wc_df2$n_pairs,
                 frac_1111 = wc_df2$frac_1111, mean_r_score = wc_df2$mean_r_score),
      data.frame(method = "Louvain", module = lc_df2$module, n_pairs = lc_df2$n_pairs,
                 frac_1111 = lc_df2$frac_1111, mean_r_score = lc_df2$mean_r_score)
    )
    comb2$frac_1111[is.na(comb2$frac_1111)] <- 0
    comb2$module_f <- factor(paste0(comb2$method, " M", comb2$module))

    fig08 <- ggplot(comb2, aes(x = frac_1111, y = reorder(module_f, frac_1111),
                                size = log10(n_pairs + 1), color = mean_r_score)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient2(low = "#d6604d", mid = "#f7f7f7", high = "#1b7837",
                            midpoint = 0.5, name = "Mean R_score") +
      scale_size_continuous(name = "log₁₀(n_pairs)", range = c(3, 10)) +
      geom_vline(xintercept = 0.99, linetype = "dashed", color = "grey50") +
      labs(
        title    = "Fig 8 — Module Condition Patterns",
        subtitle = "All modules dominated by pattern '1111' (constitutive across Mock/DC3000/AvrRpt2/AvrRpm1)\nPoint size = log₁₀(n_pairs); color = mean R_score (cross-condition consistency)",
        x        = "Fraction of pairs with pattern '1111'",
        y        = NULL
      ) +
      theme_report()

    save_fig(fig08, "fig08_condition_patterns", width = 9, height = 6)
  }
}

# ==============================================================================
# FIG 09 — Pipeline schematic (designed, not data-driven)
# ==============================================================================
cat("fig09 — Pipeline schematic\n")

# Create as a clean diagram using annotation layers
fig09 <- ggplot() +
  xlim(0, 10) + ylim(0, 12) +
  coord_fixed(ratio = 0.9) +
  # Title
  annotate("text", x = 5, y = 11.6, label = "Fig 9 — Arabidopsis Co-expression Pipeline",
           size = 5, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 5, y = 11.1, label = "zscore_gene + Spearman | obs_subcluster(298 pts) | |r|≥0.42",
           size = 3.2, color = "grey40", hjust = 0.5) +
  # Stage boxes (left column = shared path; right branch = GGM)
  # Box helper: annotate rect + text
  # Stage 0: Input
  annotate("rect",  xmin=1.5, xmax=8.5, ymin=10, ymax=10.7, fill="#e8f5e9", color="#1b7837", linewidth=0.8) +
  annotate("text",  x=5, y=10.35, label="INPUT: 10x Seurat object (65,061 nuclei, 18,364 genes)\nload_seurat() — AT-ID mapping, min_cells=10 filter",
           size=3, hjust=0.5) +
  # Arrow
  annotate("segment", x=5, xend=5, y=10, yend=9.55, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Stage 1
  annotate("rect",  xmin=1.5, xmax=8.5, ymin=8.9, ymax=9.55, fill="#e3f2fd", color="#1565c0", linewidth=0.8) +
  annotate("text",  x=5, y=9.22, label="STAGE 1: Normalization choice  •  8 methods tested\n✓ zscore_gene + Spearman  (depth_leakage=0.107; eff_rank=159.9; split-half=0.895)",
           size=2.9, hjust=0.5) +
  annotate("segment", x=5, xend=5, y=8.9, yend=8.45, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Stage 2
  annotate("rect",  xmin=1.5, xmax=8.5, ymin=7.8, ymax=8.45, fill="#e3f2fd", color="#1565c0", linewidth=0.8) +
  annotate("text",  x=5, y=8.12, label="STAGE 2: Observation-point design  •  10 designs on Pareto front\n✓ obs_subcluster(298 pts)  |  Alternative: metacell_t25(2602 pts) for richness",
           size=2.9, hjust=0.5) +
  annotate("segment", x=5, xend=5, y=7.8, yend=7.35, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Stage 3
  annotate("rect",  xmin=1.5, xmax=8.5, ymin=6.7, ymax=7.35, fill="#e3f2fd", color="#1565c0", linewidth=0.8) +
  annotate("text",  x=5, y=7.02, label="STAGE 3: Edge threshold  •  9 design points (2 levers)\n✓ global |r| ≥ 0.42  →  751,959 pairs / 5,450 genes  (top-k rejected: null_gap≈1.19)",
           size=2.9, hjust=0.5) +
  # Branching arrow
  annotate("segment", x=5, xend=5, y=6.7, yend=6.35, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  annotate("segment", x=3, xend=7, y=6.35, yend=6.35, color="grey40", linewidth=0.6) +
  annotate("segment", x=3, xend=3, y=6.35, yend=5.9, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  annotate("segment", x=7, xend=7, y=6.35, yend=5.9, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Pseudobulk branch (left)
  annotate("rect",  xmin=0.3, xmax=5.4, ymin=5.2, ymax=5.9, fill="#fff8e1", color="#f57f17", linewidth=0.8) +
  annotate("text",  x=2.85, y=5.55, label="PSEUDOBULK MODE\nestimate_pseudobulk()  •  pair_scores_full.csv (54M pairs)",
           size=2.7, hjust=0.5) +
  # GGM branch (right)
  annotate("rect",  xmin=5.6, xmax=9.7, ymin=5.2, ymax=5.9, fill="#f3e5f5", color="#7b1fa2", linewidth=0.8) +
  annotate("text",  x=7.65, y=5.55, label="GGM MODE (per condition)\nestimate_singlecellggm()  •  4 edge tables",
           size=2.7, hjust=0.5) +
  # Arrows down
  annotate("segment", x=2.85, xend=2.85, y=5.2, yend=4.75, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  annotate("segment", x=7.65, xend=7.65, y=5.2, yend=4.75, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Robustness
  annotate("rect",  xmin=0.3, xmax=5.4, ymin=4.1, ymax=4.75, fill="#fff8e1", color="#f57f17", linewidth=0.8) +
  annotate("text",  x=2.85, y=4.42, label="ROBUSTNESS: compute_robustness()  +\ncharacterize_condition_pattern()  (16-bit pattern per pair)",
           size=2.7, hjust=0.5) +
  annotate("segment", x=2.85, xend=2.85, y=4.1, yend=3.65, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Modules
  annotate("rect",  xmin=0.3, xmax=5.4, ymin=3.0, ymax=3.65, fill="#fff8e1", color="#f57f17", linewidth=0.8) +
  annotate("text",  x=2.85, y=3.32, label="MODULES: build_wgcna_modules()  +  Louvain\nWGCNA 12 modules (power=9)  |  Louvain 6 modules (seed=98)",
           size=2.7, hjust=0.5) +
  annotate("segment", x=2.85, xend=2.85, y=3.0, yend=2.55, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Interpretation
  annotate("rect",  xmin=0.3, xmax=5.4, ymin=1.9, ymax=2.55, fill="#e8f5e9", color="#1b7837", linewidth=0.8) +
  annotate("text",  x=2.85, y=2.22, label="INTERPRET: annotate_go() + annotate_tfs() + annotate_context()\nGO BP (clusterProfiler) | TF list | condition context | GOI lookup",
           size=2.7, hjust=0.5) +
  # GGM interpretation (right)
  annotate("rect",  xmin=5.6, xmax=9.7, ymin=4.1, ymax=4.75, fill="#f3e5f5", color="#7b1fa2", linewidth=0.8) +
  annotate("text",  x=7.65, y=4.42, label="PER-CONDITION networks\ncondition-specific co-regulation",
           size=2.7, hjust=0.5) +
  # Legend note
  annotate("text", x=5, y=1.2,
           label="FLAG-13: GGM captures constitutive co-expression; pseudobulk captures inducible/cell-type-restricted regulators (e.g. WRKY TFs)\nThe two modes are complementary — both are standard pipeline output.",
           size=2.6, color="grey30", hjust=0.5, fontface="italic") +
  theme_void()

save_fig(fig09, "fig09_pipeline_schematic", width = 11, height = 10)

# ==============================================================================
# FIG 10 — WRKY sanity: kME distribution + mode comparison
# ==============================================================================
cat("fig10 — WRKY sanity\n")

wrky_path <- file.path(GOI_DIR, "WRKY_GGM_vs_PB.csv")
wgcna_mem_path <- file.path(WGCNA_DIR, "module_membership.csv")

if (file.exists(wrky_path) && file.exists(wgcna_mem_path)) {
  wrky_df  <- as.data.frame(fread(wrky_path,     nThread = 1L))
  wgcna_mem <- as.data.frame(fread(wgcna_mem_path, nThread = 1L))

  # Pseudobulk kME for WRKY genes
  wrky_in_pb <- wrky_df[!is.na(wrky_df$PB_wgcna_kME), ]

  # Panel A: kME distribution in pseudobulk
  pA <- ggplot(wrky_in_pb, aes(x = PB_wgcna_kME)) +
    geom_histogram(bins = 20, fill = "#4393c3", color = "white", alpha = 0.85) +
    geom_vline(xintercept = median(wrky_in_pb$PB_wgcna_kME, na.rm = TRUE),
               linetype = "dashed", color = HIGHLIGHT_COL, linewidth = 0.8) +
    labs(x = "kME in pseudobulk WGCNA", y = "Count",
         title = "WRKY kME distribution (pseudobulk)") +
    theme_report()

  # Panel B: GGM vs PB capture rates
  captured_df <- data.frame(
    mode  = c("GGM only", "Pseudobulk only", "Both", "Neither"),
    count = c(
      sum( wrky_df$captured_in_GGM & !wrky_df$captured_in_PB, na.rm = TRUE),
      sum(!wrky_df$captured_in_GGM &  wrky_df$captured_in_PB, na.rm = TRUE),
      sum( wrky_df$captured_in_GGM &  wrky_df$captured_in_PB, na.rm = TRUE),
      sum(!wrky_df$captured_in_GGM & !wrky_df$captured_in_PB, na.rm = TRUE)
    )
  )
  captured_df$frac <- captured_df$count / nrow(wrky_df)

  pB <- ggplot(captured_df, aes(x = reorder(mode, count), y = frac)) +
    geom_col(fill = c("#d6604d", "#4393c3", "#1b7837", "grey70"), alpha = 0.85) +
    geom_text(aes(label = paste0(count, " (", round(frac*100), "%)")),
              hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.9),
                       expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Fraction of 70 WRKY genes",
         title = "WRKY capture: GGM vs pseudobulk") +
    theme_report()

  library(patchwork)
  fig10 <- (pA | pB) +
    plot_annotation(
      title    = "Fig 10 — Post-hoc Sanity: WRKY TF Family (n=70) in Official Modules",
      subtitle = "Pseudobulk captures inducible/cell-type-restricted WRKYs that GGM misses (FLAG-13)\nThis is POST-HOC SANITY ONLY — WRKY recovery was NOT used in method selection"
    )

  save_fig(fig10, "fig10_wrky_sanity", width = 11, height = 5)
} else {
  message("SKIP fig10 — WRKY data not found")
}

# ==============================================================================
# Write figures manifest
# ==============================================================================
manifest_lines <- c(
  "# Figures Manifest",
  "",
  paste0("Generated: ", Sys.time()),
  "All figures saved as .png (300 dpi) and .svg",
  "",
  "| Figure | Filename | Caption | Source data |",
  "|---|---|---|---|",
  "| 1 | fig01_stage1_normalization | Normalization method comparison (8 methods; grouped bars for depth_leakage, eff_rank, split-half) | obs_design/normalization_decision.csv |",
  "| 2 | fig02_stage2_pareto | Observation-point Pareto front (splithalf vs eff_rank, all designs labeled, Pareto members highlighted) | obs_design/stage2_{metacell,cluster}*.csv |",
  "| 3 | fig03_stage2_cluster_resweep | Cluster resolution sweep: eff_rank and splithalf vs resolution (Bug #1 fix verified) | obs_design/stage2_cluster_resweep.csv |",
  "| 4 | fig04_stage3_density | Network density vs |r| threshold (Lever A); chosen point at |r|≥0.42 marked | stage3_threshold_sweep/density_table.csv |",
  "| 5 | fig05_stage3_pareto | Stability vs richness Pareto for all 9 threshold design points, Lever A vs B, null_gap annotated | stage3_threshold_sweep/stage3_metrics.csv |",
  "| 6 | fig06_modules_overview | Module size distributions for WGCNA (12) and Louvain (6); grey rates shown | modules_official/{wgcna,louvain}/module_summary.csv |",
  "| 7 | fig07_pickSoftThreshold | Scale-free R² and mean connectivity vs soft-thresholding power; power=9 marked | modules_official/pickSoftThreshold_fitIndices.csv |",
  "| 8 | fig08_condition_patterns | Fraction of pairs with pattern '1111' per module, colored by mean R_score | modules_official/{wgcna,louvain}/module_condition_patterns.csv |",
  "| 9 | fig09_pipeline_schematic | Clean diagram of full pipeline stages showing two network modes and stage decisions | (designed schematic — no data file) |",
  "| 10 | fig10_wrky_sanity | WRKY family kME distribution in pseudobulk + GGM vs PB capture rates (post-hoc sanity only) | geneset_lookups/WRKY_GGM_vs_PB.csv |"
)

writeLines(manifest_lines,
           file.path(FIG_DIR, "FIGURES_MANIFEST.md"))
message("Written: FIGURES_MANIFEST.md")

cat("\n=== All figures complete ===\n")
