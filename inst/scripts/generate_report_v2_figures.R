#!/usr/bin/env Rscript
## generate_report_v2_figures.R
## Generates NEW figures for the v2 executive report (GGM mode + benchmark).
## Run from repo root: Rscript inst/scripts/generate_report_v2_figures.R
##
## V1 figures (pseudobulk only) are copied by the shell; this script adds:
##   fig09v2 — pipeline schematic (both modes)
##   fig11    — GGM R_score saturation
##   fig12    — GGM official module sizes (4 sets)
##   fig13    — method benchmark structural metrics
##   fig14    — cross-method ARI heatmap at thr=0.3
##   fig15    — GGM condition pattern distribution
##   fig16    — GGM vs pseudobulk comparison
##   fig17    — module condition profiles (top pathogen-activated)
##
## base-R + ggplot2 + svglite. No data.table GForce. seed=98.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(svglite)
})
if (requireNamespace("patchwork", quietly = TRUE)) library(patchwork)

set.seed(98L)

BASE_DIR    <- "results/pathogen_multiome"
FIG_DIR     <- file.path(BASE_DIR, "report/v2/figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, name, width = 9, height = 6) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot = p,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(name, ".svg")), plot = p,
         width = width, height = height, device = svglite::svglite, bg = "white")
  message("Saved: ", name)
  invisible(NULL)
}

theme_report <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#f0f0f0", color = "grey70"),
          plot.title        = element_text(face = "bold", size = 13),
          plot.subtitle     = element_text(size = 10, color = "grey40"))
}

# Colours
GGM_COL   <- "#9c27b0"   # purple = GGM mode
PB_COL    <- "#1565c0"   # blue   = pseudobulk mode
WIN_COL   <- "#1b7837"   # green  = winner/chosen
ORANGE    <- "#e65100"

cat("=== Generating v2 figures (new GGM + benchmark) ===\n\n")

# ==============================================================================
# fig09v2 — pipeline schematic (BOTH modes)
# ==============================================================================
cat("fig09v2 — pipeline schematic (both modes)\n")

fig09v2 <- ggplot() +
  xlim(0, 10) + ylim(0, 13) +
  coord_fixed(ratio = 0.9) +
  # Title
  annotate("text", x = 5, y = 12.6, label = "Arabidopsis Co-expression Pipeline — Two Complementary Modes",
           size = 4.8, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 5, y = 12.1, label = "GGM mode (per condition) | Pseudobulk mode (subcluster 298 pts)",
           size = 3.2, color = "grey40", hjust = 0.5) +
  # Shared: Input
  annotate("rect",  xmin=1, xmax=9, ymin=11, ymax=11.7, fill="#e8f5e9", color=WIN_COL, linewidth=0.8) +
  annotate("text",  x=5, y=11.35, label="INPUT: 10x Seurat object | load_seurat() | AT-ID mapping | min_cells=10",
           size=3.2, hjust=0.5) +
  annotate("segment", x=5, xend=5, y=11, yend=10.6, arrow=arrow(length=unit(0.15,"cm")), color="grey40") +
  # Shared: Stage 1 norm
  annotate("rect",  xmin=1, xmax=9, ymin=9.95, ymax=10.6, fill="#e3f2fd", color=PB_COL, linewidth=0.7) +
  annotate("text",  x=5, y=10.27, label="STAGE 1 — Normalization: zscore_gene + Spearman (depth_leakage=0.107; eff_rank=159.9)\n8 methods compared | prior-free evaluation only",
           size=2.9, hjust=0.5) +
  annotate("segment", x=3, xend=7, y=9.95, yend=9.95, color="grey40", linewidth=0.5) +
  annotate("segment", x=3, xend=3, y=9.95, yend=9.55, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  annotate("segment", x=7, xend=7, y=9.95, yend=9.55, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  # Left branch: GGM
  annotate("text", x=2.2, y=9.75, label="GGM mode", size=3.4, color=GGM_COL, fontface="bold") +
  annotate("rect",  xmin=0.2, xmax=5.5, ymin=8.7, ymax=9.55, fill="#f3e5f5", color=GGM_COL, linewidth=0.8) +
  annotate("text",  x=2.85, y=9.12, label="estimate_singlecellggm()\nPER CONDITION (4 separate networks)\n~289k–449k edges per condition | pcor≥0.02 | positive only",
           size=2.7, hjust=0.5) +
  annotate("segment", x=2.85, xend=2.85, y=8.7, yend=8.3, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  # Right branch: Pseudobulk
  annotate("text", x=7.7, y=9.75, label="Pseudobulk mode", size=3.4, color=PB_COL, fontface="bold") +
  annotate("rect",  xmin=5.6, xmax=9.8, ymin=8.7, ymax=9.55, fill="#e3f2fd", color=PB_COL, linewidth=0.8) +
  annotate("text",  x=7.7, y=9.12, label="obs_subcluster(298 pts)\nSTAGE 2 — obs-point Pareto sweep\nstability=0.895 | eff_rank=159.9",
           size=2.7, hjust=0.5) +
  annotate("segment", x=7.7, xend=7.7, y=8.7, yend=8.3, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  # GGM: robustness
  annotate("rect",  xmin=0.2, xmax=5.5, ymin=7.4, ymax=8.3, fill="#f3e5f5", color=GGM_COL, linewidth=0.8) +
  annotate("text",  x=2.85, y=7.85, label="compute_robustness() | 1,413,505 pairs\nR_score threshold: saturation at 0.3–0.5 (62,863 edges)\nR_score≥0.6 → 15,384 edges (sparse, high quality)",
           size=2.7, hjust=0.5) +
  annotate("segment", x=2.85, xend=2.85, y=7.4, yend=7.0, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  # Pseudobulk: coexpr
  annotate("rect",  xmin=5.6, xmax=9.8, ymin=7.4, ymax=8.3, fill="#e3f2fd", color=PB_COL, linewidth=0.8) +
  annotate("text",  x=7.7, y=7.85, label="estimate_pseudobulk() | 54.2M pairs\nSTAGE 3 — |r|≥0.42 chosen\nnull_gap=7.2M | 751,959 edges / 5,450 genes",
           size=2.7, hjust=0.5) +
  annotate("segment", x=7.7, xend=7.7, y=7.4, yend=7.0, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  # GGM: modules
  annotate("rect",  xmin=0.2, xmax=5.5, ymin=6.1, ymax=7.0, fill="#f3e5f5", color=GGM_COL, linewidth=0.8) +
  annotate("text",  x=2.85, y=6.55, label="METHOD BENCHMARK (6×5 grid)\n4 official module sets:\nlarge: wgcna_p1(42M/68%grey) + Louvain(14M/8%grey)\nsmall: wgcna(26M/44%grey) + Louvain(15M/2%grey)",
           size=2.5, hjust=0.5) +
  annotate("segment", x=2.85, xend=2.85, y=6.1, yend=5.7, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  # Pseudobulk: modules
  annotate("rect",  xmin=5.6, xmax=9.8, ymin=6.1, ymax=7.0, fill="#e3f2fd", color=PB_COL, linewidth=0.8) +
  annotate("text",  x=7.7, y=6.55, label="build_wgcna_modules() + Louvain\n12 WGCNA modules (power=9, grey 8.1%)\n6 Louvain modules (grey 6.6%)\nAll: constitutive ('1111') dominant pattern",
           size=2.5, hjust=0.5) +
  annotate("segment", x=7.7, xend=7.7, y=6.1, yend=5.7, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  # Shared: Interpretation
  annotate("segment", x=2.85, xend=7.7, y=5.7, yend=5.7, color="grey40", linewidth=0.5) +
  annotate("segment", x=5.275, xend=5.275, y=5.7, yend=5.3, arrow=arrow(length=unit(0.12,"cm")), color="grey40") +
  annotate("rect",  xmin=1, xmax=9, ymin=4.5, ymax=5.3, fill="#e8f5e9", color=WIN_COL, linewidth=0.8) +
  annotate("text",  x=5, y=4.9, label="INTERPRET: annotate_go() + annotate_tfs() + annotate_context(ref_condition=...)\nGO BP | TF list | condition context | GOI lookup | condition-pattern profiles",
           size=2.9, hjust=0.5) +
  # FLAG-13 note
  annotate("rect",  xmin=0.5, xmax=9.5, ymin=3.5, ymax=4.2, fill="#fff8e1", color="#f57f17", linewidth=0.6) +
  annotate("text", x=5, y=3.85, size=2.8, hjust=0.5, color="#5d3a00",
           label="FLAG-13: GGM conditions over ALL cells → misses inducible/cell-restricted TFs (e.g. most WRKY family)\nPseudobulk uses subclusters as obs → captures inducible signal. The two modes are COMPLEMENTARY.") +
  # Legend
  annotate("rect", xmin=0.5, xmax=2.8, ymin=2.5, ymax=3.1, fill="#f3e5f5", color=GGM_COL) +
  annotate("text", x=1.65, y=2.8, label="GGM mode\n(per-condition)", size=2.8, color=GGM_COL, hjust=0.5) +
  annotate("rect", xmin=3.0, xmax=5.3, ymin=2.5, ymax=3.1, fill="#e3f2fd", color=PB_COL) +
  annotate("text", x=4.15, y=2.8, label="Pseudobulk mode\n(subcluster)", size=2.8, color=PB_COL, hjust=0.5) +
  annotate("rect", xmin=5.5, xmax=7.8, ymin=2.5, ymax=3.1, fill="#e8f5e9", color=WIN_COL) +
  annotate("text", x=6.65, y=2.8, label="Shared / both\nmodes", size=2.8, color=WIN_COL, hjust=0.5) +
  annotate("rect", xmin=8.0, xmax=9.5, ymin=2.5, ymax=3.1, fill="#fff8e1", color="#f57f17") +
  annotate("text", x=8.75, y=2.8, label="Design\nnote", size=2.8, color="#5d3a00", hjust=0.5) +
  theme_void()

save_fig(fig09v2, "fig09v2_pipeline_schematic_both_modes", width = 11, height = 11)

# ==============================================================================
# fig11 — GGM R_score saturation
# ==============================================================================
cat("fig11 — GGM R_score saturation\n")

rscore_df <- data.frame(
  threshold = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
  n_edges   = c(62863, 62863, 62863, 15384, 15384, 4228),
  n_genes   = c(10358, 10358, 10358, 3441, 3441, 1651)
)
rscore_long <- rbind(
  data.frame(threshold = rscore_df$threshold, metric = "n_edges", value = rscore_df$n_edges),
  data.frame(threshold = rscore_df$threshold, metric = "n_genes", value = rscore_df$n_genes)
)

fig11 <- ggplot(rscore_df, aes(x = threshold)) +
  geom_line(aes(y = n_edges), color = GGM_COL, linewidth = 1.5) +
  geom_point(aes(y = n_edges), color = GGM_COL, size = 4) +
  # Annotate the flat region
  annotate("rect", xmin = 0.28, xmax = 0.52, ymin = -2000, ymax = 67000,
           fill = "#f3e5f5", alpha = 0.3, color = NA) +
  annotate("text", x = 0.40, y = 64000, label = "R_score 0.3, 0.4, 0.5\nIDENTICAL (62,863 edges)\nR_score is discrete: saturation",
           size = 3.2, color = GGM_COL, hjust = 0.5) +
  annotate("text", x = 0.625, y = 18000, label = "0.6: 15,384\n(same as 0.7)", size = 3, color = "grey40", hjust = 0) +
  annotate("text", x = 0.82, y = 6000, label = "0.8: 4,228", size = 3, color = "grey40", hjust = 0) +
  scale_x_continuous(breaks = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Fig 11 — GGM R_score Threshold Sweep: Network Saturation",
    subtitle = "R_score is discrete on this dataset (quartile distribution: 0%–90% = 0.25; top = 1.0)\nThresholds 0.3, 0.4, 0.5 all retain the same 62,863-edge network. Effective distinct networks = 2.",
    x        = "R_score threshold",
    y        = "Number of edges"
  ) +
  theme_report()

save_fig(fig11, "fig11_ggm_rscore_saturation", width = 8, height = 5)

# ==============================================================================
# fig12 — GGM official module sizes (4 sets)
# ==============================================================================
cat("fig12 — GGM official module sizes (4 sets)\n")

# Read the module meta files
read_mod <- function(path, set_label) {
  df <- as.data.frame(fread(path, nThread = 1L))
  df$set <- set_label
  df[, c("set", "module_id", "n_genes", "top_organ_or_condition")]
}

lw  <- tryCatch(read_mod("results/pathogen_multiome/official_modules/large_wgcna/module_meta.csv",   "large_wgcna\n(R≥0.3, WGCNA p1)"), error=function(e) NULL)
ll  <- tryCatch(read_mod("results/pathogen_multiome/official_modules/large_louvain/module_meta.csv", "large_louvain\n(R≥0.3, Louvain)"), error=function(e) NULL)
sw  <- tryCatch(read_mod("results/pathogen_multiome/official_modules/small_wgcna/module_meta.csv",   "small_wgcna\n(R≥0.6, WGCNA)"), error=function(e) NULL)
sl  <- tryCatch(read_mod("results/pathogen_multiome/official_modules/small_louvain/module_meta.csv", "small_louvain\n(R≥0.6, Louvain)"), error=function(e) NULL)

mod_df <- do.call(rbind, Filter(Negate(is.null), list(lw, ll, sw, sl)))

# Grey stats from known values
grey_stats <- data.frame(
  set       = c("large_wgcna\n(R≥0.3, WGCNA p1)", "large_louvain\n(R≥0.3, Louvain)",
                "small_wgcna\n(R≥0.6, WGCNA)", "small_louvain\n(R≥0.6, Louvain)"),
  pct_grey  = c(68.0, 8.0, 44.5, 1.9),
  n_modules = c(42, 14, 26, 15)
)

cond_colours <- c(AvrRpm1="#9c27b0", DC3000="#e65100", AvrRpt2="#1565c0", Mock="#388e3c")

if (!is.null(mod_df) && nrow(mod_df) > 0) {
  mod_df$set <- factor(mod_df$set, levels = grey_stats$set)
  mod_df$top_cond <- mod_df$top_organ_or_condition

  # Module size beeswarm/bar
  fig12a <- ggplot(mod_df, aes(x = set, y = n_genes, fill = top_cond)) +
    geom_col(data = mod_df[order(mod_df$n_genes, decreasing=TRUE), ],
             aes(group = module_id), width = 0.8, position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = cond_colours, name = "Top condition") +
    labs(x = NULL, y = "n_genes per module") +
    theme_report() +
    theme(axis.text.x = element_text(size = 8))

  # Summary metrics
  sum_df <- merge(
    aggregate(n_genes ~ set, mod_df, sum),
    grey_stats, by = "set")
  sum_df$n_modules_f <- paste0(sum_df$n_modules, " modules\n", round(sum_df$pct_grey), "% grey")
  sum_df$set <- factor(sum_df$set, levels = grey_stats$set)

  fig12 <- ggplot(sum_df, aes(x = set, y = n_genes, fill = set)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = n_modules_f), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("#ce93d8","#7b1fa2","#90caf9","#1565c0"), guide = "none") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    labs(
      title    = "Fig 12 — GGM Official Module Sets: Coverage Comparison",
      subtitle = "4 sets: 2 thresholds (R≥0.3 large / R≥0.6 small) × 2 methods (WGCNA / Louvain)\nLouvain dramatically outperforms WGCNA on gene coverage (8% vs 68% grey at R≥0.3)",
      x        = NULL,
      y        = "Genes assigned to modules"
    ) +
    theme_report() +
    theme(axis.text.x = element_text(size = 9))

  save_fig(fig12, "fig12_ggm_module_sets_overview", width = 9, height = 6)
}

# ==============================================================================
# fig13 — Method benchmark structural metrics
# ==============================================================================
cat("fig13 — method benchmark structural metrics\n")

bm_path <- "results/pathogen_multiome/method_benchmark/structural_metrics.csv"
if (file.exists(bm_path)) {
  bm <- as.data.frame(fread(bm_path, nThread = 1L))

  # Deduplicate: thresholds 0.3=0.4=0.5 and 0.6=0.7 are identical inputs
  # Keep one representative per effective group
  bm_repr <- bm[bm$threshold %in% c(0.3, 0.6), ]
  bm_repr$thr_label <- paste0("R≥", bm_repr$threshold,
                               ifelse(bm_repr$threshold == 0.3,
                                      "\n(=R≥0.4=0.5; 62,863 edges)",
                                      "\n(=R≥0.7; 15,384 edges)"))

  method_order <- c("wgcna_p1","wgcna_p4","wgcna_p6","wgcna_p8","louvain","leiden")
  method_labels <- c("WGCNA\np=1","WGCNA\np=4","WGCNA\np=6","WGCNA\np=8","Louvain","Leiden")
  bm_repr$method <- factor(bm_repr$method, levels = method_order, labels = method_labels)
  bm_repr$thr_label <- factor(bm_repr$thr_label)

  # Panel 1: pct_grey
  p_grey <- ggplot(bm_repr, aes(x = method, y = pct_grey, fill = method)) +
    geom_col(width = 0.7, alpha = 0.85) +
    facet_wrap(~ thr_label) +
    scale_fill_manual(values = c("#ce93d8","#ba68c8","#ab47bc","#8e24aa",
                                 "#1565c0","#0d47a1"), guide = "none") +
    labs(x = NULL, y = "% grey (unassigned)", title = "% Grey genes") +
    theme_report() + theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))

  # Panel 2: modularity
  p_mod <- ggplot(bm_repr, aes(x = method, y = modularity, fill = method)) +
    geom_col(width = 0.7, alpha = 0.85) +
    facet_wrap(~ thr_label) +
    scale_fill_manual(values = c("#ce93d8","#ba68c8","#ab47bc","#8e24aa",
                                 "#1565c0","#0d47a1"), guide = "none") +
    labs(x = NULL, y = "Graph modularity", title = "Graph modularity (higher = better)") +
    theme_report() + theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))

  # Panel 3: n_modules
  p_nm <- ggplot(bm_repr, aes(x = method, y = n_modules, fill = method)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_text(aes(label = n_modules), vjust = -0.3, size = 3) +
    facet_wrap(~ thr_label) +
    scale_fill_manual(values = c("#ce93d8","#ba68c8","#ab47bc","#8e24aa",
                                 "#1565c0","#0d47a1"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = NULL, y = "n modules", title = "Number of modules") +
    theme_report() + theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    fig13 <- (p_grey | p_nm | p_mod) +
      patchwork::plot_annotation(
        title    = "Fig 13 — Method Benchmark: Structural Metrics (6 methods × 2 effective thresholds)",
        subtitle = "WGCNA at R≥0.3: high grey (68–93%), few modules assigned. Louvain/Leiden: low grey (8%), all genes assigned.\nAt R≥0.6: WGCNA collapses to 3–7 modules; Louvain/Leiden remain comprehensive.\nDual-method choice: WGCNA p1 (sparse, hierarchical) + Louvain (comprehensive coverage)."
      )
  } else {
    fig13 <- p_grey
  }
  save_fig(fig13, "fig13_ggm_benchmark_structural", width = 12, height = 6)
}

# ==============================================================================
# fig14 — Cross-method ARI heatmap at thr=0.3
# ==============================================================================
cat("fig14 — cross-method ARI heatmap\n")

ari_path <- "results/pathogen_multiome/method_benchmark/ari_matrix_thr0.3.csv"
if (file.exists(ari_path)) {
  ari_mat <- as.data.frame(fread(ari_path, nThread = 1L))
  methods <- ari_mat[[1]]
  ari_mat <- ari_mat[, -1]
  rownames(ari_mat) <- methods

  # Long format
  ari_long <- data.frame(
    method1 = rep(methods, each = length(methods)),
    method2 = rep(methods, times = length(methods)),
    ARI     = as.vector(as.matrix(ari_mat))
  )
  method_labels2 <- c("WGCNA\np=1","WGCNA\np=4","WGCNA\np=6","WGCNA\np=8","Louvain","Leiden")
  ari_long$method1 <- factor(ari_long$method1, levels = methods, labels = method_labels2)
  ari_long$method2 <- factor(ari_long$method2, levels = rev(methods), labels = rev(method_labels2))

  fig14 <- ggplot(ari_long, aes(x = method1, y = method2, fill = ARI)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", ARI)), size = 3.5, color = "white") +
    scale_fill_gradient2(low = "#e3f2fd", mid = "#7b1fa2", high = "#1a237e",
                         midpoint = 0.6, limits = c(0.35, 1.0), name = "ARI") +
    labs(
      title    = "Fig 14 — Cross-Method Agreement (ARI) at R_score ≥ 0.3",
      subtitle = "Adjusted Rand Index between all method pairs. Same matrix at R≥0.4 and R≥0.5 (identical input).\nWGCNA p1 vs Louvain: 0.40 (most distinct). WGCNA p6 vs p8: 0.84 (most similar).",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"))

  save_fig(fig14, "fig14_ggm_ari_heatmap", width = 7, height = 6)
}

# ==============================================================================
# fig15 — GGM condition pattern distribution
# ==============================================================================
cat("fig15 — GGM pattern distribution\n")

pat_path <- "results/pathogen_multiome/robustness/pattern_counts.csv"
if (file.exists(pat_path)) {
  pat_df <- as.data.frame(fread(pat_path, nThread = 1L))
  pat_df$frac   <- pat_df$n_pairs / sum(pat_df$n_pairs)
  pat_df$cat    <- ifelse(grepl("^single_", pat_df$pattern_label), "single-condition",
                   ifelse(grepl("^mixed", pat_df$pattern_label), "multi-condition (mixed)",
                   ifelse(pat_df$pattern_label == "constitutive_all", "constitutive",
                   ifelse(pat_df$pattern_label %in% c("ETI_shared","pan_pathogen"), "treatment-combination",
                   "none"))))
  pat_df$label_disp <- sub("single_", "", pat_df$pattern_label)
  pat_df <- pat_df[order(pat_df$n_pairs, decreasing = TRUE), ]
  pat_df$label_disp <- factor(pat_df$label_disp, levels = rev(pat_df$label_disp))

  cat_colours <- c("single-condition" = "#7b1fa2",
                   "multi-condition (mixed)" = "#bdbdbd",
                   "treatment-combination" = "#e65100",
                   "constitutive" = WIN_COL,
                   "none" = "grey50")

  fig15 <- ggplot(pat_df, aes(x = n_pairs, y = label_disp, fill = cat)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.1f%%", frac * 100)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = cat_colours, name = "Pattern type") +
    scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    labs(
      title    = "Fig 15 — GGM Condition Pattern Distribution",
      subtitle = sprintf("1,413,505 total pairs. Single-condition pairs dominate (%.0f%%).\nConstitutive_all ('1111') = only 0.3%% — most GGM edges are condition-specific.",
                         sum(pat_df$frac[pat_df$cat == "single-condition"]) * 100),
      x        = "Number of pairs",
      y        = NULL
    ) +
    theme_report() + theme(legend.position = "bottom")

  save_fig(fig15, "fig15_ggm_condition_patterns", width = 9, height = 6)
}

# ==============================================================================
# fig16 — GGM vs pseudobulk comparison
# ==============================================================================
cat("fig16 — GGM vs pseudobulk comparison\n")

compare_df <- data.frame(
  metric    = c("Network edges", "Gene coverage", "% grey\n(official modules)",
                "Dominant\ncondition pattern",
                "Best method\ngrey rate",
                "Condition\nspecificity"),
  GGM_large  = c("62,863", "10,358", "68% (WGCNA)\n8% (Louvain)", "single-condition\n(28% AvrRpm1-only)",
                  "8%\n(Louvain)", "HIGH\n(per-condition)"),
  PB        = c("751,959", "5,450\n(at |r|≥0.42)", "8% (WGCNA)\n7% (Louvain)", "constitutive\n('1111'; 99%+ of pairs)",
                "7%\n(Louvain)", "LOW\n(single pool)")
)

# Make a styled comparison table as a plot
tab_data <- data.frame(
  row_idx = rep(seq_len(nrow(compare_df)), 2),
  col     = rep(c("GGM (R≥0.3)", "Pseudobulk (|r|≥0.42)"), each = nrow(compare_df)),
  metric  = rep(compare_df$metric, 2),
  val     = c(compare_df$GGM_large, compare_df$PB),
  mode_col = rep(c(GGM_COL, PB_COL), each = nrow(compare_df))
)

fig16 <- ggplot(tab_data, aes(x = col, y = reorder(metric, -row_idx))) +
  geom_tile(aes(fill = col), alpha = 0.12, color = "grey80") +
  geom_text(aes(label = val), size = 3.2) +
  scale_fill_manual(values = c("GGM (R≥0.3)" = GGM_COL, "Pseudobulk (|r|≥0.42)" = PB_COL),
                    guide = "none") +
  labs(
    title    = "Fig 16 — GGM vs Pseudobulk: Side-by-Side Comparison",
    subtitle = "Two complementary modes covering different aspects of co-expression.\nGGM: sparse, condition-specific. Pseudobulk: dense, constitutive co-programs.\nFLAG-13: GGM misses inducible TFs (e.g. most WRKY family); pseudobulk captures them.",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 10))

save_fig(fig16, "fig16_ggm_vs_pseudobulk", width = 9, height = 6)

# ==============================================================================
# fig17 — Module condition profiles (top pathogen-activated GGM modules)
# ==============================================================================
cat("fig17 — module condition profiles\n")

cond_path <- "results/pathogen_multiome/condition_comparison/modules_gained_in_pathogen.csv"
if (file.exists(cond_path)) {
  cond_df <- as.data.frame(fread(cond_path, nThread = 1L))

  # Take top 15 modules by max delta
  cond_df$max_delta <- apply(cond_df[, c("delta_DC3000", "delta_AvrRpt2", "delta_AvrRpm1")], 1, max)
  top15 <- head(cond_df[order(cond_df$max_delta, decreasing = TRUE), ], 15)

  top15_long <- rbind(
    data.frame(module = top15$module, condition = "Mock",    weight = top15$Mock),
    data.frame(module = top15$module, condition = "DC3000",  weight = top15$DC3000),
    data.frame(module = top15$module, condition = "AvrRpt2", weight = top15$AvrRpt2),
    data.frame(module = top15$module, condition = "AvrRpm1", weight = top15$AvrRpm1)
  )
  top15_long$condition <- factor(top15_long$condition,
                                  levels = c("Mock","DC3000","AvrRpt2","AvrRpm1"))
  top15_long$module_f <- factor(paste0("M", top15_long$module),
                                 levels = paste0("M", top15$module))

  fig17 <- ggplot(top15_long, aes(x = condition, y = module_f, fill = weight)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", weight)), size = 2.8) +
    scale_fill_gradient2(low = "#e3f2fd", mid = "#fff3e0", high = GGM_COL,
                         midpoint = median(top15_long$weight),
                         name = "Mean intra\nmodule weight") +
    labs(
      title    = "Fig 17 — GGM Module Condition Profiles (Top 15 Most Differentially Active)",
      subtitle = "Mean intramodular edge weight per condition for the 15 modules with highest treatment activation.\nModule 26 (bottom) is the only Mock-dominant module (SAR-like; delta_DC3000 = -0.023).",
      x = NULL, y = "GGM module (large_wgcna)"
    ) +
    theme_report() + theme(panel.grid = element_blank())

  save_fig(fig17, "fig17_ggm_condition_profiles", width = 8, height = 7)
}

# ==============================================================================
# Write figures manifest
# ==============================================================================
manifest <- c(
  "# Figures Manifest — v2 Report",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "## Figures carried over from v1 (pseudobulk pipeline)",
  "",
  "| Figure | Filename | Caption summary |",
  "|---|---|---|",
  "| 01 | fig01_stage1_normalization | Stage 1 normalization method comparison (8 methods) |",
  "| 02 | fig02_stage2_pareto | Stage 2 observation-point Pareto front (all designs) |",
  "| 03 | fig03_stage2_cluster_resweep | Stage 2 cluster resolution sweep (Bug #1 fixed) |",
  "| 04 | fig04_stage3_density | Stage 3 density vs |r| threshold (Lever A) |",
  "| 05 | fig05_stage3_pareto | Stage 3 stability vs richness Pareto |",
  "| 06 | fig06_modules_overview | Pseudobulk official module sizes: WGCNA(12)+Louvain(6) |",
  "| 07 | fig07_pickSoftThreshold | WGCNA pickSoftThreshold: R^2 vs power |",
  "| 08 | fig08_condition_patterns | Pseudobulk module condition patterns (all '1111') |",
  "| 10 | fig10_wrky_sanity | WRKY TF sanity: GGM vs pseudobulk capture (post-hoc) |",
  "",
  "## New figures for v2 (GGM mode + benchmark + both modes)",
  "",
  "| Figure | Filename | Caption summary | Source data |",
  "|---|---|---|---|",
  "| 09v2 | fig09v2_pipeline_schematic_both_modes | Pipeline schematic showing both GGM and pseudobulk modes | designed |",
  "| 11 | fig11_ggm_rscore_saturation | GGM R_score threshold saturation (0.3=0.4=0.5) | robustness/threshold_sweep.csv |",
  "| 12 | fig12_ggm_module_sets_overview | 4 GGM official module sets: n_modules and gene coverage | official_modules/*/module_meta.csv |",
  "| 13 | fig13_ggm_benchmark_structural | 6-method structural metrics (pct_grey, n_modules, modularity) | method_benchmark/structural_metrics.csv |",
  "| 14 | fig14_ggm_ari_heatmap | Cross-method ARI heatmap at R_score≥0.3 | method_benchmark/ari_matrix_thr0.3.csv |",
  "| 15 | fig15_ggm_condition_patterns | GGM pair condition pattern distribution | robustness/pattern_counts.csv |",
  "| 16 | fig16_ggm_vs_pseudobulk | GGM vs pseudobulk side-by-side comparison | (derived from both modes) |",
  "| 17 | fig17_ggm_condition_profiles | Top 15 pathogen-activated GGM modules (condition heatmap) | condition_comparison/modules_gained_in_pathogen.csv |"
)
writeLines(manifest, file.path(FIG_DIR, "FIGURES_MANIFEST.md"))
message("Written: FIGURES_MANIFEST.md")

cat("\n=== V2 figure generation complete ===\n")
