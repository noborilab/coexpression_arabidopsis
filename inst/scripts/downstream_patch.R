#!/usr/bin/env Rscript
# downstream_patch.R — fixes for downstream_analysis.R run
# Fixes: A1 (pfit$KS.p length-0), D1/D2 PB (column mismatch), E1 master GO,
#        Jaccard heatmap threshold, topology_degree figure

options(WGCNA.useThreads = FALSE)
suppressPackageStartupMessages({
  library(data.table); setDTthreads(1L)
  library(igraph)
  library(ggplot2)
  library(scales)
  library(base64enc)
})

REPO    <- "/Users/jep23kod/Documents/GitHub/coexpression_arabidopsis"
RESULTS <- file.path(REPO, "results/pathogen_multiome")
DOWN    <- file.path(RESULTS, "downstream")
FIGS    <- file.path(DOWN, "figures")
LOG_FILE <- file.path(REPO, "logs/downstream_overnight.log")
log_con  <- file(LOG_FILE, open="at")

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) {
  msg <- paste0("[", ts(), "] ", ..., "\n")
  cat(msg); cat(msg, file=log_con)
}
log_msg("=== downstream_patch.R START ===")

save_csv <- function(df, path) {
  write.csv(df, path, row.names=FALSE)
  log_msg("  saved ", basename(path), " (", nrow(df), " rows)")
}
save_png_svg <- function(p, stem, w=8, h=6, dpi=300) {
  ggsave(paste0(stem, ".png"), p, width=w, height=h, dpi=dpi)
  ggsave(paste0(stem, ".svg"), p, width=w, height=h)
  log_msg("  fig: ", basename(stem))
}

# ── Symbol map ────────────────────────────────────────────────────────────────
sym_map <- fread(file.path(RESULTS, "symbol_map.csv"), nThread=1L)
setnames(sym_map, c("gene_id","gene_symbol"))
sym_map <- sym_map[!duplicated(gene_id)]
label_gene <- function(ids) {
  sym <- sym_map$gene_symbol[match(ids, sym_map$gene_id)]
  ifelse(!is.na(sym) & sym != "", paste0(sym, " (", ids, ")"), ids)
}

NET_SOURCES <- list(
  Mock          = list(file=file.path(REPO,"output_per_condition/Mock/edge_table.csv"),     wt="weight"),
  DC3000        = list(file=file.path(REPO,"output_per_condition/DC3000/edge_table.csv"),   wt="weight"),
  AvrRpt2       = list(file=file.path(REPO,"output_per_condition/AvrRpt2/edge_table.csv"),  wt="weight"),
  AvrRpm1       = list(file=file.path(REPO,"output_per_condition/AvrRpm1/edge_table.csv"),  wt="weight"),
  GGM_consensus = list(file=file.path(RESULTS,"robustness/pair_scores_full.csv"),           wt="R_score", threshold=0.3),
  Pseudobulk    = list(file=file.path(RESULTS,"pseudobulk_zscore_spearman/modules_official/edges_absr042.csv"), wt="mean_abs_r")
)

# ─── FIX A1: degree distribution + power-law ──────────────────────────────────
A1_OUT <- file.path(DOWN, "topology_degree.csv")
if (!file.exists(A1_OUT)) {
  log_msg("FIX A1: degree distribution + power-law")
  rows <- list()
  for (nm in names(NET_SOURCES)) {
    src <- NET_SOURCES[[nm]]
    if (!file.exists(src$file)) next
    e <- fread(src$file, nThread=1L)
    if (!is.null(src$threshold)) e <- e[as.numeric(e[[src$wt]]) >= src$threshold, ]
    cols_ab <- grep("gene_id", names(e), value=TRUE)
    if (length(cols_ab) < 2) next
    ga <- as.character(e[[cols_ab[1]]]); gb <- as.character(e[[cols_ab[2]]])
    all_genes_vec <- c(ga, gb)
    uniq_genes    <- unique(all_genes_vec)
    deg_vec       <- tabulate(match(all_genes_vec, uniq_genes))
    pfit <- tryCatch(igraph::fit_power_law(deg_vec, xmin=NULL), error=function(e2) NULL)
    # FIX: guard against NULL or length-0 elements in pfit
    alpha_val <- if (!is.null(pfit) && !is.null(pfit$alpha)  && length(pfit$alpha)  == 1) pfit$alpha  else NA_real_
    ks_p      <- if (!is.null(pfit) && !is.null(pfit$KS.p)   && length(pfit$KS.p)   == 1) pfit$KS.p   else NA_real_
    mode_val  <- if (!is.null(src$threshold)) "ggm_consensus" else if (nm == "Pseudobulk") "pseudobulk" else "ggm_per_cond"
    rows[[nm]] <- data.frame(
      network=nm, mode=mode_val,
      gene_id=uniq_genes, degree=deg_vec,
      display_label=label_gene(uniq_genes),
      pl_alpha=alpha_val, pl_ks_p=ks_p,
      stringsAsFactors=FALSE)
    rm(e, ga, gb, all_genes_vec); gc()
    log_msg("  A1 ", nm, " done: alpha=", round(alpha_val, 3), " KS.p=", round(ks_p, 4))
  }
  deg_df <- do.call(rbind, rows)
  save_csv(deg_df, A1_OUT)
  # Figure
  p <- ggplot(deg_df, aes(x=degree, color=network)) +
    geom_density(alpha=0.6) +
    scale_x_log10(labels=comma) +
    scale_y_continuous() +
    labs(title="Degree distributions (all networks)", x="Degree (log10)", y="Density", color="Network") +
    theme_minimal(base_size=11)
  save_png_svg(p, file.path(FIGS, "fig_degree_distributions"))
  log_msg("FIX A1 done — ", nrow(deg_df), " gene-network degree records")
} else {
  log_msg("A1 already done — topology_degree.csv exists")
}

# ─── FIX D1/D2 PB: column mismatch → use rbindlist(fill=TRUE) ────────────────
PB_SETS  <- c("wgcna","louvain")
pb_mod_dir <- function(s) file.path(RESULTS, "pseudobulk_zscore_spearman/modules_official", s)

D1_PB_OUT <- file.path(DOWN, "module_condition_activation_pseudobulk.csv")
if (!file.exists(D1_PB_OUT)) {
  log_msg("FIX D1 PB: pseudobulk condition activation")
  rows_d1 <- list()
  for (s in PB_SETS) {
    mcp_f <- file.path(pb_mod_dir(s), "module_condition_patterns.csv")
    if (!file.exists(mcp_f)) next
    d <- fread(mcp_f, nThread=1L)
    d[, set := paste0("PB_",s)]
    rows_d1[[s]] <- d
  }
  if (length(rows_d1) > 0) {
    pb_act <- rbindlist(rows_d1, fill=TRUE)
    save_csv(as.data.frame(pb_act), D1_PB_OUT)
    # Figure
    pb_act2 <- as.data.frame(pb_act)
    w_cols <- c("w_Mock","w_DC3000","w_AvrRpt2","w_AvrRpm1")
    if (all(c("module","set") %in% names(pb_act2)) && any(w_cols %in% names(pb_act2))) {
      avail_w <- intersect(w_cols, names(pb_act2))
      pb_long <- reshape(pb_act2[, c("set","module", avail_w)],
                         varying=avail_w, v.names="weight", timevar="condition",
                         times=gsub("w_","",avail_w), direction="long")
      p <- ggplot(pb_long[!is.na(pb_long$weight),], aes(x=factor(module), y=weight, fill=condition)) +
        geom_col(position="dodge") + facet_wrap(~set) +
        labs(title="Pseudobulk module condition activation\n(WITHIN pseudobulk mode only; NOT cross-mode comparable)",
             x="Module", y="Mean weight") +
        theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=90, size=7))
      save_png_svg(p, file.path(FIGS, "fig_module_condition_activation_pseudobulk"), w=12, h=6)
    }
    log_msg("FIX D1 PB done — ", nrow(pb_act), " records")
  }
} else {
  log_msg("D1 PB already done")
}

D2_PB_OUT <- file.path(DOWN, "module_condition_patterns_pseudobulk.csv")
if (!file.exists(D2_PB_OUT)) {
  log_msg("FIX D2 PB: pseudobulk condition-pattern fractions")
  rows_d2 <- list()
  for (s in PB_SETS) {
    mcp_f <- file.path(pb_mod_dir(s), "module_condition_patterns.csv")
    if (!file.exists(mcp_f)) next
    d <- fread(mcp_f, nThread=1L)
    d[, set := paste0("PB_",s)]
    rows_d2[[s]] <- d
  }
  if (length(rows_d2) > 0) {
    pb_patt <- rbindlist(rows_d2, fill=TRUE)
    save_csv(as.data.frame(pb_patt), D2_PB_OUT)
    log_msg("FIX D2 PB done — ", nrow(pb_patt), " records (per-mode only)")
  }
} else {
  log_msg("D2 PB already done")
}

# ─── FIX E1: master GO file (column mismatch fix) ─────────────────────────────
MASTER_GO <- file.path(DOWN, "go_enrichment_all_sets.csv")
if (!file.exists(MASTER_GO)) {
  log_msg("FIX E1: combine GO enrichment files")
  go_files <- list.files(DOWN, pattern="^go_enrichment_", full.names=TRUE)
  go_list  <- lapply(go_files, function(f) {
    d <- tryCatch(fread(f, nThread=1L), error=function(e2) NULL)
    if (is.null(d)) return(NULL)
    # Standardize to common columns
    set_val <- gsub("go_enrichment_|\\.csv","", basename(f))
    d[, source_file := set_val]
    # Normalize column names
    if ("GO_ID" %in% names(d)) setnames(d, "GO_ID", "ID")
    if ("p.adjust" %in% names(d)) setnames(d, "p.adjust", "padj")
    keep <- intersect(c("source_file","module","ID","Description","padj"), names(d))
    d[, keep, with=FALSE]
  })
  go_list <- Filter(Negate(is.null), go_list)
  if (length(go_list) > 0) {
    all_go <- rbindlist(go_list, fill=TRUE)
    save_csv(as.data.frame(all_go), MASTER_GO)
    log_msg("FIX E1 done — ", nrow(all_go), " combined GO records")
  }
} else {
  log_msg("E1 master GO already done")
}

# ─── FIX: crossmode Jaccard heatmap (lower threshold) ────────────────────────
JACC_FIG <- file.path(FIGS, "fig_crossmode_jaccard_heatmap.png")
if (!file.exists(JACC_FIG)) {
  log_msg("FIX: crossmode Jaccard heatmap")
  ov_df <- tryCatch(read.csv(file.path(DOWN,"crossmode_overlap.csv")), error=function(e) NULL)
  if (!is.null(ov_df) && nrow(ov_df) > 0) {
    # Use all pairs or top 5% by Jaccard
    threshold <- if (max(ov_df$jaccard, na.rm=TRUE) < 0.05) 0 else 0.05
    top_ov <- ov_df[ov_df$jaccard >= threshold & ov_df$n_overlap >= 1,]
    if (nrow(top_ov) > 0) {
      # Simplify labels
      top_ov$ggm_short <- sub("GGM_[^:]+::", "GGM:", top_ov$ggm_module)
      top_ov$pb_short  <- sub("PB_[^:]+::", "PB:", top_ov$pb_module)
      top_ov$ggm_set   <- sub("::.*","", top_ov$ggm_module)
      # Show top-50 pairs by Jaccard
      top50 <- head(top_ov[order(-top_ov$jaccard),], 50)
      p <- ggplot(top50, aes(x=ggm_short, y=pb_short, fill=jaccard, size=n_overlap)) +
        geom_point(shape=21) +
        scale_fill_gradient(low="white", high="steelblue") +
        scale_size_continuous(range=c(1,8)) +
        labs(title=paste0("GGM vs Pseudobulk module overlap (top 50 pairs by Jaccard)\n",
                          "Max Jaccard=", round(max(ov_df$jaccard),4),
                          "; low overlap expected (different edge definitions)"),
             x="GGM module", y="Pseudobulk module", fill="Jaccard", size="n_overlap") +
        theme_minimal(base_size=8) +
        theme(axis.text.x=element_text(angle=90, hjust=1, size=6),
              axis.text.y=element_text(size=6))
      save_png_svg(p, file.path(FIGS, "fig_crossmode_jaccard_heatmap"), w=12, h=8)
    }
  }
  log_msg("FIX Jaccard heatmap done")
}

# ─── Regenerate FIGURES_MANIFEST.md ──────────────────────────────────────────
fig_files <- list.files(FIGS, pattern="\\.(png|svg)$", full.names=FALSE)
fig_manifest <- c(
  "# Figures Manifest",
  paste0("Generated: ", ts()),
  "",
  "| Filename | Caption | Source CSV |",
  "|----------|---------|------------|"
)
captions <- list(
  "fig_degree_distributions.png"       = list(c="Degree distributions (log10) for all 6 networks; power-law fit alpha per network", s="topology_degree.csv"),
  "fig_hub_genes.png"                  = list(c="Centrality scatter: degree vs eigenvector (NA for large networks, degree-only)", s="topology_centrality.csv"),
  "fig_global_stats_comparison.png"    = list(c="Global network statistics per network (GGM and pseudobulk)", s="topology_global_stats.csv"),
  "fig_module_kme.png"                 = list(c="kME distributions per module set (all 6 sets)", s="module_kme_distributions.csv"),
  "fig_module_quality_across_sets.png" = list(c="Module quality: grey rate, n_modules, kME median across all 6 sets", s="module_quality_summary.csv"),
  "fig_crossmode_jaccard_heatmap.png"  = list(c="GGM vs pseudobulk module overlap (Jaccard); low overlap expected due to different edge definitions", s="crossmode_overlap.csv"),
  "fig_crossmethod_ari.png"            = list(c="WGCNA vs Louvain ARI within each mode (GGM large, GGM small, pseudobulk)", s="crossmethod_agreement.csv"),
  "fig_core_modespecific.png"          = list(c="Core (both modes) vs GGM-only vs pseudobulk-only gene counts", s="core_vs_modespecific.csv"),
  "fig_cross_set_consistency.png"      = list(c="GGM cross-set stability: histogram of genes by number of GGM sets assigned", s="cross_set_consistency.csv"),
  "fig_module_condition_activation_ggm.png"        = list(c="GGM module condition activation (within GGM only; not cross-mode comparable)", s="module_condition_activation_ggm.csv"),
  "fig_module_condition_activation_pseudobulk.png" = list(c="Pseudobulk module condition activation (within PB only; not cross-mode comparable)", s="module_condition_activation_pseudobulk.csv"),
  "fig_wrky_crossmode.png"             = list(c="WRKY family enrichment by module — demo of generic gene-set query", s="geneset_query_demo_wrky.csv")
)
for (f in sort(fig_files)) {
  info <- captions[[f]]
  cap  <- if (!is.null(info)) info$c else "Supplementary figure"
  src  <- if (!is.null(info)) info$s else "—"
  fig_manifest <- c(fig_manifest, paste0("| ", f, " | ", cap, " | ", src, " |"))
}
writeLines(fig_manifest, file.path(DOWN, "FIGURES_MANIFEST.md"))
log_msg("FIGURES_MANIFEST.md updated — ", length(fig_files), " figures")

# ─── Final status ─────────────────────────────────────────────────────────────
outputs <- c("topology_degree.csv","topology_global_stats.csv","topology_centrality.csv",
             "module_kme_distributions.csv","module_eigengene_correlations.csv",
             "module_hubs.csv","module_quality_summary.csv",
             "crossmode_overlap.csv","crossmode_ari.csv","crossmethod_agreement.csv",
             "core_vs_modespecific.csv","cross_set_consistency.csv",
             "module_condition_activation_ggm.csv","module_condition_activation_pseudobulk.csv",
             "module_condition_patterns_ggm.csv","module_condition_patterns_pseudobulk.csv",
             "go_enrichment_all_sets.csv","tf_enrichment.csv",
             "gene_lookup_master.csv","geneset_query_demo_wrky.csv",
             "DOWNSTREAM_ANALYSIS_REPORT.html")
for (f in outputs) {
  exists_str <- if (file.exists(file.path(DOWN,f))) "DONE" else "MISSING"
  log_msg("  ", exists_str, ": ", f)
}
log_msg("=== downstream_patch.R COMPLETE ===")
close(log_con)
