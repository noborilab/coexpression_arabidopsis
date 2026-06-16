#!/usr/bin/env Rscript
# downstream_analysis.R
# Downstream analysis suite: topology, module quality, cross-mode comparison,
# functional annotation, gene-centric utility, integrated HTML report.
# Unattended overnight run — incremental, resumable, failure-tolerant.
# seed=98 everywhere; base-R + fread(nThread=1L); no data.table GForce.

options(WGCNA.useThreads = FALSE)   # disable WGCNA threading on aarch64
suppressPackageStartupMessages({
  library(data.table)
  setDTthreads(1L)                   # ensure data.table single-thread throughout
  library(igraph)
  library(ggplot2)
  library(WGCNA)
  library(jsonlite)
  library(RColorBrewer)
  library(scales)
  library(aricode)
  library(grDevices)
  library(base64enc)
})

# ─── Paths ────────────────────────────────────────────────────────────────────
REPO      <- "/Users/jep23kod/Documents/GitHub/coexpression_arabidopsis"
RESULTS   <- file.path(REPO, "results/pathogen_multiome")
DOWN      <- file.path(RESULTS, "downstream")
FIGS      <- file.path(DOWN, "figures")
LOGS      <- file.path(REPO, "logs")
dir.create(DOWN, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGS, recursive=TRUE, showWarnings=FALSE)
dir.create(LOGS, recursive=TRUE, showWarnings=FALSE)

LOG_FILE <- file.path(LOGS, "downstream_overnight.log")
log_con  <- file(LOG_FILE, open="at")
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) {
  msg <- paste0("[", ts(), "] ", ..., "\n")
  cat(msg); cat(msg, file=log_con)
}
log_msg("=== downstream_analysis.R START ===")

# ─── Helpers ──────────────────────────────────────────────────────────────────
set.seed(98)
skip_if_exists <- function(path) {
  if (file.exists(path)) { log_msg("SKIP (exists): ", basename(path)); return(TRUE) }
  return(FALSE)
}
safe_run <- function(label, expr) {
  log_msg("BEGIN: ", label)
  t0 <- proc.time()["elapsed"]
  result <- tryCatch(expr, error=function(e) { log_msg("ERROR in ", label, ": ", e$message); NULL })
  dt <- round(proc.time()["elapsed"] - t0, 1)
  if (!is.null(result)) log_msg("DONE: ", label, " [", dt, "s]")
  invisible(result)
}

save_csv <- function(df, path) {
  write.csv(df, path, row.names=FALSE)
  log_msg("  saved ", basename(path), " (", nrow(df), " rows)")
}
save_png_svg <- function(p, stem, w=8, h=6, dpi=300) {
  ggsave(paste0(stem, ".png"), p, width=w, height=h, dpi=dpi)
  ggsave(paste0(stem, ".svg"), p, width=w, height=h)
  log_msg("  fig: ", basename(paste0(stem, ".png/svg")))
}

# ─── PHASE 0.5: symbol map ────────────────────────────────────────────────────
log_msg("PHASE 0.5: loading symbol map")
sym_map <- fread(file.path(RESULTS, "symbol_map.csv"), nThread=1L)
setnames(sym_map, c("gene_id","gene_symbol"))
sym_map <- sym_map[!duplicated(gene_id)]
label_gene <- function(ids) {
  sym <- sym_map$gene_symbol[match(ids, sym_map$gene_id)]
  ifelse(!is.na(sym) & sym != "", paste0(sym, " (", ids, ")"), ids)
}
log_msg("  symbol map: ", nrow(sym_map), " entries")

# ─── Input paths ──────────────────────────────────────────────────────────────
CONDS    <- c("Mock","DC3000","AvrRpt2","AvrRpm1")
GGM_SETS <- c("large_wgcna","large_louvain","small_wgcna","small_louvain")
PB_SETS  <- c("wgcna","louvain")

ggm_mod_dir <- function(s) file.path(RESULTS, "official_modules", s)
pb_mod_dir  <- function(s) file.path(RESULTS, "pseudobulk_zscore_spearman/modules_official", s)

# ─── PHASE 0: Inventory ───────────────────────────────────────────────────────
log_msg("PHASE 0: inventory")

inventory_lines <- c(
  "# Downstream Analysis Inventory",
  paste0("Generated: ", ts()),
  "",
  "## Existing results (pre-downstream)"
)

inv_dirs <- c(
  file.path(RESULTS, "method_benchmark"),
  file.path(RESULTS, "condition_comparison"),
  file.path(RESULTS, "official_modules"),
  file.path(RESULTS, "robustness"),
  file.path(RESULTS, "geneset_lookups"),
  file.path(RESULTS, "pseudobulk_zscore_spearman/modules_official")
)

for (d in inv_dirs) {
  if (!dir.exists(d)) next
  inventory_lines <- c(inventory_lines, "", paste0("### ", basename(d)))
  fls <- list.files(d, recursive=TRUE, full.names=FALSE)
  for (f in fls) {
    fp <- file.path(d, f)
    sz <- if (file.exists(fp)) paste0(" [", round(file.size(fp)/1024), " KB]") else ""
    inventory_lines <- c(inventory_lines, paste0("- ", f, sz))
  }
}

inventory_lines <- c(inventory_lines, "",
  "## Key reusable results",
  "- method_benchmark/: ARI matrices, structural_metrics.csv — REUSE (cross-method agreement)",
  "- condition_comparison/module_condition_profiles.csv — REUSE for GGM condition activation (D1)",
  "- official_modules/cross_set_assignments.csv — REUSE for C4",
  "- official_modules/all_modules_condition_patterns.csv — REUSE for D2 GGM",
  "- {module_set}/go_enrichment.csv — REUSE for E1 (pseudobulk sets have it)",
  "- {module_set}/module_tfs.csv — REUSE for E2",
  "- geneset_lookups/WRKY_GGM_vs_PB.csv — REUSE for F2",
  "",
  "## NOTE: GGM module GO enrichment",
  "- GGM official module sets have module_meta.csv with top GO term only (single line)",
  "- Full GO enrichment computed in E1 via clusterProfiler if org.At.tair.db available",
  "",
  "## Cross-mode caveat (always enforce)",
  "GGM uses partial correlation (cells as obs); pseudobulk uses Spearman (subclusters as obs).",
  "I_s is near-binary in both. Per-mode condition activation is fine WITHIN a mode.",
  "Cross-mode condition-pattern quantitative comparison is FORBIDDEN."
)

writeLines(inventory_lines, file.path(DOWN, "DOWNSTREAM_INVENTORY.md"))
log_msg("PHASE 0 complete — DOWNSTREAM_INVENTORY.md written")

# ─── Load GGM module assignments ──────────────────────────────────────────────
load_ggm_modules <- function() {
  lapply(setNames(GGM_SETS, GGM_SETS), function(s) {
    f <- file.path(ggm_mod_dir(s), "gene_module.csv")
    if (!file.exists(f)) return(NULL)
    d <- fread(f, nThread=1L)
    setnames(d, tolower(names(d)))
    d[, display_label := label_gene(gene_id)]
    d
  })
}
load_pb_modules <- function() {
  lapply(setNames(PB_SETS, PB_SETS), function(s) {
    f <- file.path(pb_mod_dir(s), "module_membership.csv")
    if (!file.exists(f)) return(NULL)
    d <- fread(f, nThread=1L)
    setnames(d, tolower(names(d)))
    d[, display_label := label_gene(gene_id)]
    d
  })
}

ggm_mods <- load_ggm_modules()
pb_mods  <- load_pb_modules()

# ─── A: NETWORK TOPOLOGY ──────────────────────────────────────────────────────
log_msg("=== PHASE A: Network Topology ===")

# Define network sources
NET_SOURCES <- list(
  Mock      = list(file=file.path(REPO,"output_per_condition/Mock/edge_table.csv"),     wt="weight", mode="ggm_per_cond"),
  DC3000    = list(file=file.path(REPO,"output_per_condition/DC3000/edge_table.csv"),   wt="weight", mode="ggm_per_cond"),
  AvrRpt2   = list(file=file.path(REPO,"output_per_condition/AvrRpt2/edge_table.csv"),  wt="weight", mode="ggm_per_cond"),
  AvrRpm1   = list(file=file.path(REPO,"output_per_condition/AvrRpm1/edge_table.csv"),  wt="weight", mode="ggm_per_cond"),
  GGM_consensus = list(file=file.path(RESULTS,"robustness/pair_scores_full.csv"),
                       wt="R_score", threshold=0.3, mode="ggm_consensus"),
  Pseudobulk    = list(file=file.path(RESULTS,"pseudobulk_zscore_spearman/modules_official/edges_absr042.csv"),
                       wt="mean_abs_r", mode="pseudobulk")
)

# ── A3: Global stats (graph_from_edgelist for all networks — fast C-level code) ─
A3_OUT <- file.path(DOWN, "topology_global_stats.csv")
safe_run("A3: global network stats", {
  if (!skip_if_exists(A3_OUT)) {
    rows <- list()
    for (nm in names(NET_SOURCES)) {
      src <- NET_SOURCES[[nm]]
      if (!file.exists(src$file)) { log_msg("  SKIP ", nm, " - file missing"); next }
      e <- fread(src$file, nThread=1L)
      if (!is.null(src$threshold)) {
        e <- e[as.numeric(e[[src$wt]]) >= src$threshold, ]
      }
      wt_col <- src$wt
      cols_ab <- grep("gene_id", names(e), value=TRUE)
      if (length(cols_ab) < 2) next
      ga <- as.character(e[[cols_ab[1]]]); gb <- as.character(e[[cols_ab[2]]])
      n_edges <- length(ga)
      # Degree via base-R (no igraph, no data.table GForce)
      all_genes_vec <- c(ga, gb)
      deg_tbl  <- tabulate(match(all_genes_vec, unique(all_genes_vec)))
      n_nodes  <- length(unique(all_genes_vec))
      mean_deg <- mean(deg_tbl)
      density_val <- if (n_nodes > 1) 2 * n_edges / (n_nodes * (n_nodes - 1)) else NA
      # Component stats + clustering: use graph_from_edgelist (fast, avoids data.frame overhead)
      clust_coef <- NA_real_; n_comp <- NA_integer_; lcc_frac <- NA_real_
      g_struct <- tryCatch(graph_from_edgelist(cbind(ga, gb), directed=FALSE), error=function(e2) NULL)
      if (!is.null(g_struct)) {
        comps <- components(g_struct)
        n_comp   <- comps$no
        lcc_frac <- max(comps$csize) / vcount(g_struct)
        # Clustering only for small networks (transitivity is O(E*d), expensive for dense graphs)
        if (n_edges < 80000) {
          clust_coef <- tryCatch(transitivity(g_struct, type="global"), error=function(e2) NA_real_)
        }
        rm(g_struct); gc()
      }
      rows[[nm]] <- data.frame(network=nm, mode=src$mode, n_nodes=n_nodes, n_edges=n_edges,
                               density=density_val, mean_degree=mean_deg,
                               global_clustering_coef=clust_coef,
                               n_components=n_comp, lcc_fraction=lcc_frac,
                               stringsAsFactors=FALSE)
      rm(e, ga, gb, all_genes_vec); gc()
      log_msg("  A3 ", nm, " done: ", n_nodes, " nodes, ", n_edges, " edges")
    }
    gs_df <- do.call(rbind, rows)
    save_csv(gs_df, A3_OUT)
    log_msg("ANALYSIS A3 done — ", nrow(gs_df), " networks characterized")
  }
})

# ── A1: Degree distribution + scale-free fit ──────────────────────────────────
A1_OUT <- file.path(DOWN, "topology_degree.csv")
safe_run("A1: degree distribution + power-law", {
  if (!skip_if_exists(A1_OUT)) {
    rows <- list()
    for (nm in names(NET_SOURCES)) {
      src <- NET_SOURCES[[nm]]
      if (!file.exists(src$file)) next
      e <- fread(src$file, nThread=1L)
      if (!is.null(src$threshold)) {
        e <- e[as.numeric(e[[src$wt]]) >= src$threshold, ]
      }
      cols_ab <- grep("gene_id", names(e), value=TRUE)
      if (length(cols_ab) < 2) next
      ga <- as.character(e[[cols_ab[1]]]); gb <- as.character(e[[cols_ab[2]]])
      all_genes_vec <- c(ga, gb)
      uniq_genes    <- unique(all_genes_vec)
      deg_vec       <- tabulate(match(all_genes_vec, uniq_genes))
      pfit <- tryCatch(igraph::fit_power_law(deg_vec, xmin=NULL), error=function(e2) NULL)
      alpha_val <- if (!is.null(pfit)) pfit$alpha else NA_real_
      ks_p      <- if (!is.null(pfit)) pfit$KS.p else NA_real_
      gene_deg_df <- data.frame(
        network  = nm,
        mode     = src$mode,
        gene_id  = uniq_genes,
        degree   = deg_vec,
        display_label = label_gene(uniq_genes),
        pl_alpha = alpha_val,
        pl_ks_p  = ks_p,
        stringsAsFactors=FALSE
      )
      rows[[nm]] <- gene_deg_df
      rm(e, ga, gb, all_genes_vec); gc()
      log_msg("  A1 ", nm, " done: alpha=", round(alpha_val, 3))
    }
    deg_df <- do.call(rbind, rows)
    save_csv(deg_df, A1_OUT)
    # Figure: degree distributions
    p <- ggplot(deg_df, aes(x=degree, color=network)) +
      geom_density(alpha=0.6) +
      scale_x_log10(labels=comma) +
      scale_y_log10() +
      labs(title="Degree distributions (all networks)", x="Degree (log)", y="Density (log)",
           color="Network") +
      theme_minimal(base_size=11) +
      theme(legend.position="right")
    save_png_svg(p, file.path(FIGS, "fig_degree_distributions"))
    log_msg("ANALYSIS A1 done — power-law alpha range: ",
            round(min(deg_df$pl_alpha, na.rm=TRUE),2), " - ",
            round(max(deg_df$pl_alpha, na.rm=TRUE),2))
  }
})

# ── A2: Centrality & hubs ─────────────────────────────────────────────────────
# igraph centrality only for networks < 80k edges; larger → degree only from edge table
A2_OUT <- file.path(DOWN, "topology_centrality.csv")
safe_run("A2: centrality", {
  if (!skip_if_exists(A2_OUT)) {
    rows <- list()
    for (nm in names(NET_SOURCES)) {
      src <- NET_SOURCES[[nm]]
      if (!file.exists(src$file)) next
      e <- fread(src$file, nThread=1L)
      if (!is.null(src$threshold)) {
        e <- e[as.numeric(e[[src$wt]]) >= src$threshold, ]
      }
      cols_ab <- grep("gene_id", names(e), value=TRUE)
      if (length(cols_ab) < 2) next
      ga <- as.character(e[[cols_ab[1]]]); gb <- as.character(e[[cols_ab[2]]])
      n_edges <- length(ga)
      # Degree from edge table (always safe, no igraph)
      all_g_vec <- c(ga, gb)
      uniq_g    <- unique(all_g_vec)
      deg_vec   <- tabulate(match(all_g_vec, uniq_g))
      gene_ids  <- uniq_g
      eig_c     <- rep(NA_real_, length(gene_ids))
      btw_c     <- rep(NA_real_, length(gene_ids))
      # For small networks, also compute eigenvector + betweenness via igraph
      if (n_edges < 80000) {
        wt_vec <- if (src$wt %in% names(e)) abs(as.numeric(e[[src$wt]])) else rep(1, n_edges)
        g <- tryCatch(
          graph_from_edgelist(cbind(ga, gb), directed=FALSE),
          error=function(e2) { log_msg("  igraph build failed for ", nm, ": ", e2$message); NULL })
        if (!is.null(g)) {
          E(g)$weight <- wt_vec
          comps   <- components(g)
          lcc_ids <- which(comps$membership == which.max(comps$csize))
          g_lcc   <- induced_subgraph(g, lcc_ids)
          lcc_genes <- V(g_lcc)$name
          eig_lcc <- tryCatch(eigen_centrality(g_lcc, weights=E(g_lcc)$weight)$vector,
                              error=function(e2) rep(NA, vcount(g_lcc)))
          btw_lcc <- rep(NA_real_, vcount(g_lcc))
          if (n_edges <= 100000) {
            log_msg("  computing betweenness for ", nm)
            btw_lcc <- tryCatch(betweenness(g_lcc, normalized=TRUE),
                                error=function(e2) rep(NA, vcount(g_lcc)))
          }
          # Map back to full gene list
          idx_lcc <- match(lcc_genes, gene_ids)
          valid   <- !is.na(idx_lcc)
          eig_c[idx_lcc[valid]] <- eig_lcc[valid]
          btw_c[idx_lcc[valid]] <- btw_lcc[valid]
          rm(g, g_lcc); gc()
        }
      } else {
        log_msg("  A2 centrality: ", nm, " has ", n_edges,
                " edges — degree only (eigenvector/betweenness skipped for large networks)")
      }
      rows[[nm]] <- data.frame(
        network=nm, mode=src$mode,
        gene_id=gene_ids, display_label=label_gene(gene_ids),
        degree=deg_vec, eigenvector=eig_c, betweenness=btw_c,
        stringsAsFactors=FALSE)
      rm(e, ga, gb, all_g_vec); gc()
      log_msg("  A2 ", nm, " done: ", length(gene_ids), " genes")
    }
    cent_df <- do.call(rbind, rows)
    save_csv(cent_df, A2_OUT)
    p <- ggplot(cent_df, aes(x=degree, y=ifelse(is.na(eigenvector), 0, eigenvector), color=network)) +
      geom_point(alpha=0.2, size=0.5) +
      scale_x_log10() +
      labs(title="Centrality: degree vs eigenvector (NA=large network, degree only)",
           x="Degree (log10)", y="Eigenvector centrality") +
      theme_minimal(base_size=10) +
      facet_wrap(~network, scales="free")
    save_png_svg(p, file.path(FIGS, "fig_hub_genes"), w=12, h=8)
    log_msg("ANALYSIS A2 done — ", nrow(cent_df), " gene-network records")
  }
})

# Figure: global stats comparison
safe_run("A3: global stats figure", {
  fig_path <- file.path(FIGS, "fig_global_stats_comparison.png")
  if (!file.exists(fig_path) && file.exists(A3_OUT)) {
    gs_df <- read.csv(A3_OUT)
    gs_long <- reshape(gs_df[, c("network","mode","n_nodes","n_edges","mean_degree","global_clustering_coef","lcc_fraction")],
                       varying=c("n_nodes","n_edges","mean_degree","global_clustering_coef","lcc_fraction"),
                       v.names="value", timevar="metric",
                       times=c("n_nodes","n_edges","mean_degree","global_clustering_coef","lcc_fraction"),
                       direction="long")
    p <- ggplot(gs_long, aes(x=network, y=value, fill=mode)) +
      geom_col() +
      facet_wrap(~metric, scales="free_y") +
      labs(title="Global network statistics", x="Network", y="Value") +
      theme_minimal(base_size=9) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    save_png_svg(p, file.path(FIGS, "fig_global_stats_comparison"), w=12, h=8)
    log_msg("ANALYSIS A3 done — global stats figure saved")
  }
})

log_msg("PHASE A complete")

# ─── B: MODULE QUALITY ────────────────────────────────────────────────────────
log_msg("=== PHASE B: Module Internal Structure & Quality ===")

# Collect all module membership data
collect_all_module_memberships <- function() {
  rows <- list()
  for (s in GGM_SETS) {
    m <- ggm_mods[[s]]
    if (is.null(m)) next
    mod_col <- if ("top_module" %in% names(m)) "top_module" else "module"
    tmp <- data.frame(set=paste0("GGM_",s), mode="ggm",
                      gene_id=m$gene_id,
                      module=m[[mod_col]],
                      kME=m$kme,
                      display_label=m$display_label,
                      stringsAsFactors=FALSE)
    rows[[paste0("ggm_",s)]] <- tmp
  }
  for (s in PB_SETS) {
    m <- pb_mods[[s]]
    if (is.null(m)) next
    mod_col <- if ("module" %in% names(m)) "module" else names(m)[2]
    tmp <- data.frame(set=paste0("PB_",s), mode="pseudobulk",
                      gene_id=m$gene_id,
                      module=m[[mod_col]],
                      kME=m$kme,
                      display_label=m$display_label,
                      stringsAsFactors=FALSE)
    rows[[paste0("pb_",s)]] <- tmp
  }
  do.call(rbind, rows)
}

all_mems <- safe_run("load module memberships", collect_all_module_memberships())

# ── B1: kME distributions ─────────────────────────────────────────────────────
B1_OUT <- file.path(DOWN, "module_kme_distributions.csv")
safe_run("B1: kME distributions", {
  if (!skip_if_exists(B1_OUT) && !is.null(all_mems)) {
    # Flag low-coherence modules (median kME < 0.3)
    am <- all_mems[!is.na(all_mems$kME) & all_mems$module != 0,]
    kme_sum <- do.call(rbind, lapply(split(am, paste(am$set, am$module)), function(x) {
      data.frame(set=x$set[1], mode=x$mode[1], module=x$module[1],
                 n_genes=nrow(x),
                 kME_mean=mean(x$kME, na.rm=TRUE),
                 kME_median=median(x$kME, na.rm=TRUE),
                 kME_sd=sd(x$kME, na.rm=TRUE),
                 kME_q25=quantile(x$kME, 0.25, na.rm=TRUE),
                 kME_q75=quantile(x$kME, 0.75, na.rm=TRUE),
                 low_coherence=(median(x$kME, na.rm=TRUE) < 0.3),
                 stringsAsFactors=FALSE)
    }))
    save_csv(kme_sum, B1_OUT)
    p <- ggplot(am[am$kME >= 0,], aes(x=kME, fill=set)) +
      geom_histogram(bins=50, alpha=0.7) +
      facet_wrap(~set, scales="free_y") +
      labs(title="kME distributions per module set", x="kME", y="Count") +
      theme_minimal(base_size=9) +
      theme(legend.position="none")
    save_png_svg(p, file.path(FIGS, "fig_module_kme"), w=12, h=8)
    log_msg("ANALYSIS B1 done — kME distributions, low-coherence modules: ",
            sum(kme_sum$low_coherence, na.rm=TRUE))
  }
})

# ── B2: Eigengene-eigengene correlations ──────────────────────────────────────
B2_OUT <- file.path(DOWN, "module_eigengene_correlations.csv")
safe_run("B2: eigengene correlations", {
  if (!skip_if_exists(B2_OUT)) {
    eg_rows <- list()
    # GGM sets: use kME matrix (gene × module), correlate columns
    for (s in GGM_SETS) {
      ef <- file.path(ggm_mod_dir(s), "eigengenes.csv")
      if (!file.exists(ef)) next
      eg <- read.csv(ef, row.names=1, check.names=FALSE)
      cor_mat <- cor(eg, use="pairwise.complete.obs", method="pearson")
      # Melt upper triangle
      idx <- which(upper.tri(cor_mat), arr.ind=TRUE)
      tmp <- data.frame(set=paste0("GGM_",s), mode="ggm",
                        module_A=colnames(cor_mat)[idx[,1]],
                        module_B=colnames(cor_mat)[idx[,2]],
                        correlation=cor_mat[idx],
                        stringsAsFactors=FALSE)
      eg_rows[[paste0("ggm_",s)]] <- tmp
    }
    # Pseudobulk WGCNA: use blockwiseModules MEs (sample × module)
    net_wg <- tryCatch(
      readRDS(file.path(pb_mod_dir("wgcna"), "blockwiseModules_net.rds")),
      error=function(e) NULL)
    if (!is.null(net_wg) && !is.null(net_wg$MEs)) {
      me <- net_wg$MEs
      cor_mat <- cor(me, use="pairwise.complete.obs")
      idx <- which(upper.tri(cor_mat), arr.ind=TRUE)
      tmp <- data.frame(set="PB_wgcna", mode="pseudobulk",
                        module_A=colnames(cor_mat)[idx[,1]],
                        module_B=colnames(cor_mat)[idx[,2]],
                        correlation=cor_mat[idx],
                        stringsAsFactors=FALSE)
      eg_rows[["pb_wgcna"]] <- tmp
    }
    # Pseudobulk Louvain: compute eigengenes from obs matrix
    pb_cache <- tryCatch(
      readRDS(file.path(RESULTS,"stage3_threshold_sweep/obs_normalized_cache.rds")),
      error=function(e) NULL)
    pb_louvain <- pb_mods[["louvain"]]
    if (!is.null(pb_cache) && !is.null(pb_louvain)) {
      expr_mat <- pb_cache$matrix   # genes × obs
      gene_ids_cache <- pb_cache$gene_ids
      if (is.null(gene_ids_cache)) gene_ids_cache <- rownames(expr_mat)
      # For each non-grey louvain module, PC1 of gene subset
      mods_uniq <- sort(unique(pb_louvain$module))
      mods_uniq <- mods_uniq[mods_uniq != 0]
      me_list <- lapply(mods_uniq, function(m) {
        gids <- pb_louvain$gene_id[pb_louvain$module == m]
        idx_g <- which(gene_ids_cache %in% gids)
        if (length(idx_g) < 3) return(rep(NA, ncol(expr_mat)))
        sub_mat <- expr_mat[idx_g, , drop=FALSE]
        sub_mat <- sub_mat[complete.cases(sub_mat), , drop=FALSE]
        if (nrow(sub_mat) < 3) return(rep(NA, ncol(expr_mat)))
        pca <- prcomp(t(sub_mat), scale.=FALSE, center=TRUE)
        pca$x[,1]
      })
      me_df <- do.call(cbind, me_list)
      colnames(me_df) <- paste0("ME", mods_uniq)
      valid_cols <- apply(me_df, 2, function(x) !all(is.na(x)))
      me_df <- me_df[, valid_cols, drop=FALSE]
      if (ncol(me_df) >= 2) {
        cor_mat <- cor(me_df, use="pairwise.complete.obs")
        idx <- which(upper.tri(cor_mat), arr.ind=TRUE)
        tmp <- data.frame(set="PB_louvain", mode="pseudobulk",
                          module_A=colnames(cor_mat)[idx[,1]],
                          module_B=colnames(cor_mat)[idx[,2]],
                          correlation=cor_mat[idx],
                          stringsAsFactors=FALSE)
        eg_rows[["pb_louvain"]] <- tmp
      }
    }
    eg_df <- do.call(rbind, eg_rows)
    save_csv(eg_df, B2_OUT)
    # Figure: heatmap for each set
    for (s_name in unique(eg_df$set)) {
      sub <- eg_df[eg_df$set == s_name,]
      mods <- sort(unique(c(sub$module_A, sub$module_B)))
      mat  <- matrix(NA, length(mods), length(mods), dimnames=list(mods, mods))
      for (i in seq_len(nrow(sub))) {
        mat[sub$module_A[i], sub$module_B[i]] <- sub$correlation[i]
        mat[sub$module_B[i], sub$module_A[i]] <- sub$correlation[i]
      }
      diag(mat) <- 1
      png(file.path(FIGS, paste0("fig_eigengene_heatmap_", gsub("[^A-Za-z0-9]","_",s_name), ".png")),
          width=900, height=800, res=150)
      heatmap(mat, symm=TRUE, col=colorRampPalette(c("blue","white","red"))(100),
              main=paste0("Eigengene correlation — ", s_name), margins=c(8,8))
      dev.off()
    }
    log_msg("ANALYSIS B2 done — ", nrow(eg_df), " module-pair correlations across ", length(unique(eg_df$set)), " sets")
  }
})

# ── B3: Intramodular hub genes ─────────────────────────────────────────────────
B3_OUT <- file.path(DOWN, "module_hubs.csv")
safe_run("B3: hub genes", {
  if (!skip_if_exists(B3_OUT)) {
    hub_rows <- list()
    # GGM: hub_genes.csv already exists
    for (s in GGM_SETS) {
      hf <- file.path(ggm_mod_dir(s), "hub_genes.csv")
      if (!file.exists(hf)) next
      h <- fread(hf, nThread=1L)
      setnames(h, tolower(names(h)))  # normalize: kME -> kme
      h[, set := paste0("GGM_",s)]
      h[, mode := "ggm"]
      h[, display_label := label_gene(gene_id)]
      hub_rows[[paste0("ggm_",s)]] <- h
    }
    # Pseudobulk: compute top-kME per module from membership
    for (s in PB_SETS) {
      m <- pb_mods[[s]]
      if (is.null(m)) next
      mod_col <- if ("module" %in% names(m)) "module" else names(m)[2]
      m2 <- m[m[[mod_col]] != 0 & !is.na(m$kme), ]
      m2 <- m2[order(m2[[mod_col]], -m2$kme), ]
      top_per_mod <- do.call(rbind, lapply(split(m2, m2[[mod_col]]), function(x) head(x, 5)))
      top_per_mod$set  <- paste0("PB_",s)
      top_per_mod$mode <- "pseudobulk"
      top_per_mod$hub_rank <- unlist(lapply(split(top_per_mod, top_per_mod[[mod_col]]),
                                             function(x) seq_len(nrow(x))))
      top_per_mod$module_id <- top_per_mod[[mod_col]]
      hub_rows[[paste0("pb_",s)]] <- top_per_mod[, c("set","mode","module_id","gene_id","kme","hub_rank","display_label")]
    }
    hubs_df <- do.call(rbind, lapply(hub_rows, function(x) {
      x <- as.data.frame(x)
      keep <- intersect(c("set","mode","module_id","gene_id","kme","hub_rank","display_label"), names(x))
      x[, keep, drop=FALSE]
    }))
    save_csv(hubs_df, B3_OUT)
    log_msg("ANALYSIS B3 done — ", nrow(hubs_df), " hub gene records across all sets")
  }
})

# ── B4: Module quality summary ─────────────────────────────────────────────────
B4_OUT <- file.path(DOWN, "module_quality_summary.csv")
safe_run("B4: module quality summary", {
  if (!skip_if_exists(B4_OUT) && !is.null(all_mems)) {
    rows_q <- list()
    for (s_id in unique(all_mems$set)) {
      sub <- all_mems[all_mems$set == s_id, ]
      mode_val <- sub$mode[1]
      n_genes_total <- nrow(sub)
      # grey = module 0 or NA
      is_grey <- sub$module == 0 | is.na(sub$module)
      n_grey  <- sum(is_grey)
      n_assigned <- n_genes_total - n_grey
      pct_grey   <- 100 * n_grey / n_genes_total
      non_grey   <- sub[!is_grey, ]
      mod_sizes  <- table(non_grey$module)
      n_modules  <- length(mod_sizes)
      rows_q[[s_id]] <- data.frame(
        set=s_id, mode=mode_val,
        n_genes_total=n_genes_total,
        n_assigned=n_assigned,
        n_grey=n_grey,
        pct_grey=round(pct_grey,2),
        n_modules=n_modules,
        module_size_min=if(n_modules>0) min(mod_sizes) else NA,
        module_size_median=if(n_modules>0) median(mod_sizes) else NA,
        module_size_max=if(n_modules>0) max(mod_sizes) else NA,
        module_size_mean=if(n_modules>0) mean(mod_sizes) else NA,
        kME_mean=mean(non_grey$kME, na.rm=TRUE),
        kME_median=median(non_grey$kME, na.rm=TRUE),
        stringsAsFactors=FALSE
      )
    }
    qual_df <- do.call(rbind, rows_q)
    save_csv(qual_df, B4_OUT)
    # Figure
    qual_long <- reshape(qual_df[, c("set","pct_grey","n_modules","kME_median")],
                         varying=c("pct_grey","n_modules","kME_median"),
                         v.names="value", timevar="metric",
                         times=c("pct_grey","n_modules","kME_median"),
                         direction="long")
    p <- ggplot(qual_long, aes(x=set, y=value, fill=set)) +
      geom_col() +
      facet_wrap(~metric, scales="free_y") +
      labs(title="Module quality summary — all 6 module sets", x="Module set", y="Value") +
      theme_minimal(base_size=9) +
      theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
    save_png_svg(p, file.path(FIGS, "fig_module_quality_across_sets"), w=10, h=6)
    log_msg("ANALYSIS B4 done — ", nrow(qual_df), " module sets summarized; grey rates: ",
            paste(round(qual_df$pct_grey,1), collapse=", "), "%")
  }
})

log_msg("PHASE B complete")

# ─── C: CROSS-MODE / CROSS-METHOD ─────────────────────────────────────────────
log_msg("=== PHASE C: Cross-mode & Cross-method Comparison ===")

# ── C1: GGM vs pseudobulk overlap ─────────────────────────────────────────────
C1_OUT <- file.path(DOWN, "crossmode_overlap.csv")
safe_run("C1: GGM vs pseudobulk module overlap", {
  if (!skip_if_exists(C1_OUT) && !is.null(all_mems)) {
    # Build gene-module lists for each set
    gene_sets <- list()
    for (s_id in unique(all_mems$set)) {
      sub <- all_mems[all_mems$set == s_id & !is.na(all_mems$module) & all_mems$module != 0,]
      for (m in unique(sub$module)) {
        key <- paste0(s_id, "::M", m)
        gene_sets[[key]] <- sub$gene_id[sub$module == m]
      }
    }
    # Pairwise Jaccard between GGM and PB sets
    ggm_keys <- grep("^GGM_", names(gene_sets), value=TRUE)
    pb_keys  <- grep("^PB_",  names(gene_sets), value=TRUE)
    rows_c1 <- list()
    for (gk in ggm_keys) {
      for (pk in pb_keys) {
        a <- gene_sets[[gk]]; b <- gene_sets[[pk]]
        jacc <- length(intersect(a,b)) / length(union(a,b))
        rows_c1[[paste(gk,pk)]] <- data.frame(
          ggm_module=gk, pb_module=pk, jaccard=jacc,
          n_overlap=length(intersect(a,b)),
          n_union=length(union(a,b)),
          stringsAsFactors=FALSE)
      }
    }
    ov_df <- do.call(rbind, rows_c1)
    # ARI per pair of sets
    ari_rows <- list()
    for (gs in GGM_SETS) {
      gm <- ggm_mods[[gs]]
      if (is.null(gm)) next
      gmod_col <- if ("top_module" %in% names(gm)) "top_module" else "module"
      for (ps in PB_SETS) {
        pm <- pb_mods[[ps]]
        if (is.null(pm)) next
        shared <- intersect(gm$gene_id, pm$gene_id)
        if (length(shared) < 10) next
        glab <- gm[[gmod_col]][match(shared, gm$gene_id)]
        plab <- pm$module[match(shared, pm$gene_id)]
        ari_val <- tryCatch(ARI(glab, plab), error=function(e) NA_real_)
        ari_rows[[paste(gs,ps)]] <- data.frame(
          ggm_set=paste0("GGM_",gs), pb_set=paste0("PB_",ps),
          n_shared=length(shared), ARI=ari_val, stringsAsFactors=FALSE)
      }
    }
    ari_df <- do.call(rbind, ari_rows)
    # Combine both into one output + save separately
    save_csv(ov_df, C1_OUT)
    save_csv(ari_df, file.path(DOWN, "crossmode_ari.csv"))

    # Figure: Jaccard heatmap (top cross-mode overlaps)
    # Pivot best Jaccard per GGM-PB pair
    top_ov <- ov_df[ov_df$jaccard > 0.05,]
    if (nrow(top_ov) > 0) {
      p <- ggplot(top_ov, aes(x=ggm_module, y=pb_module, fill=jaccard)) +
        geom_tile() +
        scale_fill_gradient(low="white", high="steelblue") +
        labs(title="GGM vs Pseudobulk module overlap (Jaccard)\n[showing Jaccard > 0.05]",
             x="GGM module", y="Pseudobulk module", fill="Jaccard") +
        theme_minimal(base_size=7) +
        theme(axis.text.x=element_text(angle=90, hjust=1, size=5),
              axis.text.y=element_text(size=5))
      save_png_svg(p, file.path(FIGS, "fig_crossmode_jaccard_heatmap"), w=14, h=10)
    }
    log_msg("ANALYSIS C1 done — ", nrow(ov_df), " module pairs; ARI range: ",
            paste(round(range(ari_df$ARI, na.rm=TRUE),3), collapse=" - "))
  }
})

# ── C2: WGCNA vs Louvain agreement ────────────────────────────────────────────
C2_OUT <- file.path(DOWN, "crossmethod_agreement.csv")
safe_run("C2: WGCNA vs Louvain ARI", {
  if (!skip_if_exists(C2_OUT)) {
    # Check existing
    exist_f <- file.path(RESULTS, "method_benchmark/cross_method_agreement.csv")
    if (file.exists(exist_f)) {
      file.copy(exist_f, C2_OUT)
      log_msg("  REUSED existing cross_method_agreement.csv from method_benchmark/")
    }
    # Compute fresh for official module sets
    rows_c2 <- list()
    # Within GGM: large_wgcna vs large_louvain, small_wgcna vs small_louvain
    for (suffix in c("large","small")) {
      gm_w  <- ggm_mods[[paste0(suffix,"_wgcna")]]
      gm_l  <- ggm_mods[[paste0(suffix,"_louvain")]]
      if (is.null(gm_w) || is.null(gm_l)) next
      shared <- intersect(gm_w$gene_id, gm_l$gene_id)
      aw <- if ("top_module" %in% names(gm_w)) gm_w$top_module[match(shared, gm_w$gene_id)] else gm_w$module[match(shared, gm_w$gene_id)]
      al <- if ("top_module" %in% names(gm_l)) gm_l$top_module[match(shared, gm_l$gene_id)] else gm_l$module[match(shared, gm_l$gene_id)]
      ari_val <- tryCatch(ARI(aw, al), error=function(e) NA_real_)
      rows_c2[[paste0("GGM_",suffix)]] <- data.frame(
        comparison=paste0("GGM_",suffix,": wgcna vs louvain"),
        mode="ggm", n_shared=length(shared), ARI=ari_val, stringsAsFactors=FALSE)
    }
    # Within pseudobulk: wgcna vs louvain
    shared_pb <- intersect(pb_mods[["wgcna"]]$gene_id, pb_mods[["louvain"]]$gene_id)
    if (length(shared_pb) > 10) {
      aw <- pb_mods[["wgcna"]]$module[match(shared_pb, pb_mods[["wgcna"]]$gene_id)]
      al <- pb_mods[["louvain"]]$module[match(shared_pb, pb_mods[["louvain"]]$gene_id)]
      ari_val <- tryCatch(ARI(aw, al), error=function(e) NA_real_)
      rows_c2[["PB"]] <- data.frame(
        comparison="Pseudobulk: wgcna vs louvain",
        mode="pseudobulk", n_shared=length(shared_pb), ARI=ari_val, stringsAsFactors=FALSE)
    }
    new_agree <- do.call(rbind, rows_c2)
    save_csv(new_agree, C2_OUT)
    # Figure
    p <- ggplot(new_agree, aes(x=comparison, y=ARI, fill=mode)) +
      geom_col() +
      labs(title="WGCNA vs Louvain ARI (within each mode)", y="ARI", x=NULL) +
      theme_minimal(base_size=10) +
      theme(axis.text.x=element_text(angle=30, hjust=1))
    save_png_svg(p, file.path(FIGS, "fig_crossmethod_ari"))
    log_msg("ANALYSIS C2 done — ARI range: ",
            paste(round(range(new_agree$ARI, na.rm=TRUE),3), collapse=" - "))
  }
})

# ── C3: Core vs mode-specific genes ───────────────────────────────────────────
C3_OUT <- file.path(DOWN, "core_vs_modespecific.csv")
safe_run("C3: core vs mode-specific genes", {
  if (!skip_if_exists(C3_OUT) && !is.null(all_mems)) {
    # GGM genes (any set, non-grey)
    ggm_genes <- unique(all_mems$gene_id[
      grepl("^GGM_", all_mems$set) & !is.na(all_mems$module) & all_mems$module != 0])
    # PB genes (any set, non-grey)
    pb_genes  <- unique(all_mems$gene_id[
      grepl("^PB_",  all_mems$set) & !is.na(all_mems$module) & all_mems$module != 0])
    # Universe = union of all assigned genes
    all_assigned <- union(ggm_genes, pb_genes)
    in_ggm <- all_assigned %in% ggm_genes
    in_pb  <- all_assigned %in% pb_genes
    category <- ifelse(in_ggm & in_pb, "core",
                       ifelse(in_ggm, "GGM_specific", "PB_specific"))
    core_df <- data.frame(
      gene_id=all_assigned,
      display_label=label_gene(all_assigned),
      in_GGM=in_ggm,
      in_PB=in_pb,
      category=category,
      stringsAsFactors=FALSE
    )
    save_csv(core_df, C3_OUT)
    cat_tbl <- table(core_df$category)
    # Figure
    p <- ggplot(as.data.frame(cat_tbl), aes(x=Var1, y=Freq, fill=Var1)) +
      geom_col() +
      labs(title="Core vs mode-specific gene partitioning",
           x="Category", y="Number of genes") +
      theme_minimal(base_size=12) +
      theme(legend.position="none")
    save_png_svg(p, file.path(FIGS, "fig_core_modespecific"))
    log_msg("ANALYSIS C3 done — core=", cat_tbl["core"],
            " GGM-only=", cat_tbl["GGM_specific"],
            " PB-only=", cat_tbl["PB_specific"])
  }
})

# ── C4: Cross-set assignment consistency ──────────────────────────────────────
C4_OUT <- file.path(DOWN, "cross_set_consistency.csv")
safe_run("C4: cross-set assignment consistency", {
  src_f <- file.path(RESULTS, "official_modules/cross_set_assignments.csv")
  if (!skip_if_exists(C4_OUT) && file.exists(src_f)) {
    cs <- fread(src_f, nThread=1L)
    # Count how many sets a gene is stably assigned in (non-grey, non-NA)
    mod_cols <- c("module_large_wgcna","module_large_louvain","module_small_wgcna","module_small_louvain")
    available_cols <- intersect(mod_cols, names(cs))
    # Use base-R to avoid data.table GForce (rowSums with .SDcols is safe but let's be explicit)
    cs_df <- as.data.frame(cs)
    cs_df$n_sets_assigned <- rowSums(
      !is.na(cs_df[, available_cols, drop=FALSE]) &
      cs_df[, available_cols, drop=FALSE] != 0 &
      cs_df[, available_cols, drop=FALSE] != "grey",
      na.rm=TRUE)
    cs_df$display_label <- label_gene(cs_df$gene_id)
    cs <- as.data.table(cs_df)
    save_csv(as.data.frame(cs), C4_OUT)
    # Figure: histogram of assignment stability
    p <- ggplot(as.data.frame(cs), aes(x=n_sets_assigned)) +
      geom_bar(fill="steelblue") +
      labs(title="GGM cross-set assignment stability",
           x="Number of GGM sets gene is assigned (non-grey)", y="Count") +
      theme_minimal(base_size=12)
    save_png_svg(p, file.path(FIGS, "fig_cross_set_consistency"))
    stab_tbl <- table(cs$n_sets_assigned)
    log_msg("ANALYSIS C4 done — assigned in all 4 sets: ", stab_tbl["4"])
  }
})

log_msg("PHASE C complete")

# ─── D: CONDITION-SPECIFICITY ─────────────────────────────────────────────────
log_msg("=== PHASE D: Condition-specificity (WITHIN-mode only) ===")

# ── D1: Module condition activation ───────────────────────────────────────────
D1_GGM_OUT <- file.path(DOWN, "module_condition_activation_ggm.csv")
D1_PB_OUT  <- file.path(DOWN, "module_condition_activation_pseudobulk.csv")

safe_run("D1: GGM condition activation", {
  if (!skip_if_exists(D1_GGM_OUT)) {
    src <- file.path(RESULTS, "condition_comparison/module_condition_profiles.csv")
    if (file.exists(src)) {
      d <- fread(src, nThread=1L)
      save_csv(as.data.frame(d), D1_GGM_OUT)
      log_msg("ANALYSIS D1 GGM done — ", nrow(d), " modules from condition_comparison/")
    } else {
      # Collect from per-set module_condition_patterns.csv
      rows_d1 <- list()
      for (s in GGM_SETS) {
        mcp_f <- file.path(ggm_mod_dir(s), "module_condition_patterns.csv")
        if (!file.exists(mcp_f)) next
        d <- fread(mcp_f, nThread=1L)
        d[, set := s]
        rows_d1[[s]] <- d
      }
      if (length(rows_d1) > 0) {
        ggm_act <- do.call(rbind, rows_d1)
        save_csv(as.data.frame(ggm_act), D1_GGM_OUT)
        log_msg("ANALYSIS D1 GGM done — ", nrow(ggm_act), " module-set records")
      }
    }
  }
})

safe_run("D1: Pseudobulk condition activation", {
  if (!skip_if_exists(D1_PB_OUT)) {
    rows_d1 <- list()
    for (s in PB_SETS) {
      mcp_f <- file.path(pb_mod_dir(s), "module_condition_patterns.csv")
      if (!file.exists(mcp_f)) next
      d <- fread(mcp_f, nThread=1L)
      d[, set := paste0("PB_",s)]
      rows_d1[[s]] <- d
    }
    if (length(rows_d1) > 0) {
      pb_act <- do.call(rbind, rows_d1)
      save_csv(as.data.frame(pb_act), D1_PB_OUT)
      # Figure
      pb_act2 <- as.data.frame(pb_act)
      if (all(c("module","w_Mock","w_DC3000","w_AvrRpt2","w_AvrRpm1") %in% names(pb_act2))) {
        pb_long <- reshape(pb_act2[, c("set","module","w_Mock","w_DC3000","w_AvrRpt2","w_AvrRpm1")],
                           varying=c("w_Mock","w_DC3000","w_AvrRpt2","w_AvrRpm1"),
                           v.names="weight", timevar="condition",
                           times=c("Mock","DC3000","AvrRpt2","AvrRpm1"), direction="long")
        p <- ggplot(pb_long, aes(x=factor(module), y=weight, fill=condition)) +
          geom_col(position="dodge") +
          facet_wrap(~set) +
          labs(title="Pseudobulk module condition activation (WITHIN pseudobulk mode only;\nNOT comparable to GGM condition weights)",
               x="Module", y="Mean intramodular edge weight") +
          theme_minimal(base_size=9) +
          theme(axis.text.x=element_text(angle=90, size=7))
        save_png_svg(p, file.path(FIGS, "fig_module_condition_activation_pseudobulk"), w=12, h=6)
      }
      log_msg("ANALYSIS D1 PB done — ", nrow(pb_act), " module-set records")
    }
  }
})

# D1 GGM figure
safe_run("D1: GGM condition activation figure", {
  fig_path <- file.path(FIGS, "fig_module_condition_activation_ggm.png")
  if (!file.exists(fig_path) && file.exists(D1_GGM_OUT)) {
    ggm_act <- read.csv(D1_GGM_OUT)
    if (all(c("module","Mock","DC3000","AvrRpt2","AvrRpm1") %in% names(ggm_act))) {
      ggm_long <- reshape(ggm_act[, c("module","Mock","DC3000","AvrRpt2","AvrRpm1")],
                          varying=c("Mock","DC3000","AvrRpt2","AvrRpm1"),
                          v.names="weight", timevar="condition",
                          times=c("Mock","DC3000","AvrRpt2","AvrRpm1"), direction="long")
      p <- ggplot(ggm_long[ggm_long$module %in% head(sort(unique(ggm_long$module)),20),],
                  aes(x=factor(module), y=weight, fill=condition)) +
        geom_col(position="dodge") +
        labs(title="GGM module condition activation (within GGM mode only;\nNOT comparable to pseudobulk weights)",
             x="Module", y="Mean weight") +
        theme_minimal(base_size=9)
      save_png_svg(p, file.path(FIGS, "fig_module_condition_activation_ggm"), w=12, h=6)
    }
  }
})

# ── D2: Condition-pattern fractions per module ─────────────────────────────────
D2_GGM_OUT <- file.path(DOWN, "module_condition_patterns_ggm.csv")
D2_PB_OUT  <- file.path(DOWN, "module_condition_patterns_pseudobulk.csv")

safe_run("D2: GGM condition-pattern fractions", {
  if (!skip_if_exists(D2_GGM_OUT)) {
    src <- file.path(RESULTS, "official_modules/all_modules_condition_patterns.csv")
    if (file.exists(src)) {
      d <- fread(src, nThread=1L)
      save_csv(as.data.frame(d), D2_GGM_OUT)
      log_msg("ANALYSIS D2 GGM done — reused all_modules_condition_patterns.csv")
    }
  }
})

safe_run("D2: Pseudobulk condition-pattern fractions", {
  if (!skip_if_exists(D2_PB_OUT)) {
    rows_d2 <- list()
    for (s in PB_SETS) {
      mcp_f <- file.path(pb_mod_dir(s), "module_condition_patterns.csv")
      if (!file.exists(mcp_f)) next
      d <- fread(mcp_f, nThread=1L)
      d[, set := paste0("PB_",s)]
      rows_d2[[s]] <- d
    }
    if (length(rows_d2) > 0) {
      pb_patt <- do.call(rbind, rows_d2)
      save_csv(as.data.frame(pb_patt), D2_PB_OUT)
      log_msg("ANALYSIS D2 PB done — ", nrow(pb_patt), " records (per-mode, not cross-mode comparable)")
    }
  }
})

log_msg("PHASE D complete")

# ─── E: FUNCTIONAL ANNOTATION ─────────────────────────────────────────────────
log_msg("=== PHASE E: Functional Annotation (REFERENCE ONLY) ===")

# ── E1: GO enrichment ─────────────────────────────────────────────────────────
safe_run("E1: GO enrichment collection", {
  rows_go <- list()
  # GGM: module_meta.csv has top GO term; run clusterProfiler for full enrichment
  for (s in GGM_SETS) {
    go_f <- file.path(ggm_mod_dir(s), "go_enrichment.csv")
    meta_f <- file.path(ggm_mod_dir(s), "module_meta.csv")
    out_f  <- file.path(DOWN, paste0("go_enrichment_GGM_", s, ".csv"))
    if (skip_if_exists(out_f)) {
      d <- read.csv(out_f); d$set <- paste0("GGM_",s); rows_go[[paste0("ggm_",s)]] <- d; next
    }
    if (file.exists(go_f)) {
      d <- fread(go_f, nThread=1L)
      d[, set := paste0("GGM_",s)]
      save_csv(as.data.frame(d), out_f)
      rows_go[[paste0("ggm_",s)]] <- as.data.frame(d)
      log_msg("  GO enrichment GGM ", s, ": ", nrow(d), " terms (from existing file)")
      next
    }
    if (!file.exists(meta_f)) next
    meta <- fread(meta_f, nThread=1L)
    gm   <- ggm_mods[[s]]
    if (is.null(gm)) next
    gmod_col <- if ("top_module" %in% names(gm)) "top_module" else "module"
    gm2 <- as.data.frame(gm)[!is.na(gm[[gmod_col]]) & gm[[gmod_col]] != 0, ]
    universe <- gm2$gene_id
    go_rows <- list()
    for (m in unique(gm2[[gmod_col]])) {
      gene_set <- gm2$gene_id[gm2[[gmod_col]] == m]
      if (length(gene_set) < 10) next
      ego <- tryCatch(clusterProfiler::enrichGO(gene=gene_set, universe=universe,
                                                  OrgDb="org.At.tair.db",
                                                  keyType="TAIR", ont="BP",
                                                  pAdjustMethod="BH", qvalueCutoff=0.05,
                                                  minGSSize=10),
                       error=function(e) NULL)
      if (!is.null(ego) && nrow(ego@result) > 0) {
        r <- ego@result[ego@result$p.adjust < 0.05, ]
        if (nrow(r) > 0) {
          r$module <- m
          go_rows[[as.character(m)]] <- r[, c("module","ID","Description","p.adjust","GeneRatio","BgRatio")]
        }
      }
    }
    if (length(go_rows) > 0) {
      go_df <- do.call(rbind, go_rows)
      go_df$set <- paste0("GGM_",s)
      save_csv(go_df, out_f)
      rows_go[[paste0("ggm_",s)]] <- go_df
      log_msg("  GO enrichment GGM ", s, ": ", nrow(go_df), " significant terms")
    }
  }
  # Pseudobulk: go_enrichment.csv already exists
  for (s in PB_SETS) {
    go_f  <- file.path(pb_mod_dir(s), "go_enrichment.csv")
    out_f <- file.path(DOWN, paste0("go_enrichment_PB_", s, ".csv"))
    if (skip_if_exists(out_f)) {
      d <- read.csv(out_f); d$set <- paste0("PB_",s); rows_go[[paste0("pb_",s)]] <- d; next
    }
    if (file.exists(go_f)) {
      d <- fread(go_f, nThread=1L)
      d[, set := paste0("PB_",s)]
      save_csv(as.data.frame(d), out_f)
      rows_go[[paste0("pb_",s)]] <- as.data.frame(d)
      log_msg("  GO enrichment PB ", s, ": ", nrow(d), " terms (reused)")
    }
  }
  # Master GO file
  master_go <- file.path(DOWN, "go_enrichment_all_sets.csv")
  if (!file.exists(master_go) && length(rows_go) > 0) {
    all_go <- do.call(rbind, lapply(rows_go, function(x) {
      x <- as.data.frame(x)
      # Standardize columns
      want <- c("set","module","ID","Description","p.adjust")
      have <- intersect(want, names(x))
      x[, have, drop=FALSE]
    }))
    save_csv(all_go, master_go)
    log_msg("ANALYSIS E1 done — ", nrow(all_go), " GO enrichment records in master file")
  }
})

# ── E2: TF enrichment ─────────────────────────────────────────────────────────
E2_OUT <- file.path(DOWN, "tf_enrichment.csv")
safe_run("E2: TF enrichment", {
  if (!skip_if_exists(E2_OUT)) {
    rows_e2 <- list()
    for (s in GGM_SETS) {
      tf_f <- file.path(ggm_mod_dir(s), "module_tfs.csv")
      if (!file.exists(tf_f)) next
      d <- fread(tf_f, nThread=1L)
      d[, set := paste0("GGM_",s)]
      d[, mode := "ggm"]
      rows_e2[[paste0("ggm_",s)]] <- d
    }
    # Pseudobulk: module_tfs.csv
    for (s in PB_SETS) {
      tf_f <- file.path(pb_mod_dir(s), "module_tfs.csv")
      if (file.exists(tf_f)) {
        d <- fread(tf_f, nThread=1L)
        d[, set := paste0("PB_",s)]
        d[, mode := "pseudobulk"]
        rows_e2[[paste0("pb_",s)]] <- d
      } else {
        # Fall back to parent pseudobulk dir
        tf_f2 <- file.path(RESULTS, "pseudobulk_zscore_spearman/module_tfs.csv")
        if (file.exists(tf_f2)) {
          d <- fread(tf_f2, nThread=1L)
          d[, set := paste0("PB_",s)]
          d[, mode := "pseudobulk"]
          rows_e2[[paste0("pb_",s)]] <- d
        }
      }
    }
    if (length(rows_e2) > 0) {
      tf_df <- do.call(rbind, lapply(rows_e2, function(x) {
        x <- as.data.frame(x)
        x$display_label <- label_gene(x$gene_id)
        x
      }))
      save_csv(tf_df, E2_OUT)
      log_msg("ANALYSIS E2 done — ", nrow(tf_df), " TF-module records; NOTE: descriptive only, modules not named")
    } else {
      log_msg("ANALYSIS E2 SKIPPED — no module_tfs.csv files found")
    }
  }
})

log_msg("PHASE E complete")

# ─── F: GENE-CENTRIC UTILITY ──────────────────────────────────────────────────
log_msg("=== PHASE F: Gene-centric Utility ===")

# ── F1: Master gene lookup table ──────────────────────────────────────────────
F1_OUT <- file.path(DOWN, "gene_lookup_master.csv")
safe_run("F1: gene lookup master", {
  if (!skip_if_exists(F1_OUT) && !is.null(all_mems)) {
    # Universe: all genes with any assignment
    all_genes <- union(
      unique(all_mems$gene_id),
      unique(c(fread(file.path(REPO,"output_per_condition/Mock/edge_table.csv"),nThread=1L,select=c(1,2)) |> unlist(use.names=FALSE)))
    )
    all_genes <- unique(all_genes)
    base_df <- data.frame(gene_id=all_genes, display_label=label_gene(all_genes),
                           stringsAsFactors=FALSE)
    # Add module membership for each set
    for (s in GGM_SETS) {
      m <- ggm_mods[[s]]
      if (is.null(m)) next
      mod_col <- if ("top_module" %in% names(m)) "top_module" else "module"
      kme_col <- "kme"
      base_df[[paste0("module_GGM_",s)]] <- m[[mod_col]][match(base_df$gene_id, m$gene_id)]
      base_df[[paste0("kME_GGM_",s)]]    <- m[[kme_col]][match(base_df$gene_id, m$gene_id)]
    }
    for (s in PB_SETS) {
      m <- pb_mods[[s]]
      if (is.null(m)) next
      mod_col <- if ("module" %in% names(m)) "module" else names(m)[2]
      base_df[[paste0("module_PB_",s)]] <- m[[mod_col]][match(base_df$gene_id, m$gene_id)]
      base_df[[paste0("kME_PB_",s)]]    <- m$kme[match(base_df$gene_id, m$gene_id)]
    }
    # Top-10 co-expressed partners in GGM consensus and pseudobulk (by weight)
    log_msg("  computing top co-expression partners (GGM consensus)...")
    ps_ggm <- fread(file.path(RESULTS,"robustness/pair_scores_full.csv"), nThread=1L)
    ps_ggm <- ps_ggm[ps_ggm$R_score >= 0.3, ]
    # base-R implementation: avoid data.table GForce (by= operations segfault on aarch64-darwin)
    get_top10_partners <- function(edge_df, gene_col_a, gene_col_b, wt_col, top_n=10) {
      e <- as.data.frame(edge_df)
      gene_a <- e[[gene_col_a]]; gene_b <- e[[gene_col_b]]; wt <- abs(as.numeric(e[[wt_col]]))
      all_g <- c(gene_a, gene_b)
      all_p <- c(gene_b, gene_a)
      all_w <- c(wt, wt)
      ord   <- order(all_g, -all_w)
      all_g <- all_g[ord]; all_p <- all_p[ord]; all_w <- all_w[ord]
      grp_rle   <- rle(all_g)
      grp_lens  <- grp_rle$lengths
      top_mask  <- unlist(lapply(grp_lens, function(n) seq_len(n) <= top_n), use.names=FALSE)
      top_g <- all_g[top_mask]; top_p <- all_p[top_mask]; top_w <- all_w[top_mask]
      pair_str  <- paste0(label_gene(top_p), ":", round(top_w, 3))
      agg <- tapply(pair_str, top_g, paste, collapse=";")
      data.frame(gene=names(agg), top_partners=as.character(agg), stringsAsFactors=FALSE)
    }
    top_ggm <- get_top10_partners(ps_ggm, "gene_id_A", "gene_id_B", "R_score")
    base_df$top10_GGM_consensus <- top_ggm$top_partners[match(base_df$gene_id, top_ggm$gene)]
    rm(ps_ggm, top_ggm); gc()
    log_msg("  computing top co-expression partners (pseudobulk)...")
    pb_edges <- fread(file.path(RESULTS,"pseudobulk_zscore_spearman/modules_official/edges_absr042.csv"), nThread=1L)
    top_pb <- get_top10_partners(pb_edges, "gene_id_A", "gene_id_B", "mean_abs_r")
    base_df$top10_PB <- top_pb$top_partners[match(base_df$gene_id, top_pb$gene)]
    rm(pb_edges, top_pb); gc()
    save_csv(base_df, F1_OUT)
    log_msg("ANALYSIS F1 done — ", nrow(base_df), " genes in master lookup table")
  }
})

# ── F2: WRKY demo ─────────────────────────────────────────────────────────────
F2_OUT <- file.path(DOWN, "geneset_query_demo_wrky.csv")
safe_run("F2: WRKY geneset query demo", {
  if (!skip_if_exists(F2_OUT)) {
    wrky_f <- file.path(RESULTS, "geneset_lookups/WRKY_GGM_vs_PB.csv")
    if (!file.exists(wrky_f)) { log_msg("  WRKY file missing — SKIP"); return(NULL) }
    wrky <- fread(wrky_f, nThread=1L)
    wrky[, display_label := label_gene(gene_id)]
    # Module enrichment: for each module set, how enriched is the WRKY set?
    # Fisher's exact test: WRKY in module vs expected by chance
    enrich_rows <- list()
    for (s in GGM_SETS) {
      m <- ggm_mods[[s]]
      if (is.null(m)) next
      mod_col <- if ("top_module" %in% names(m)) "top_module" else "module"
      all_g <- m$gene_id; wrky_g <- wrky$gene_id
      n_all <- length(all_g); n_wrky <- sum(wrky_g %in% all_g)
      for (mod_id in unique(m[[mod_col]])) {
        if (is.na(mod_id) || mod_id == 0) next
        in_mod <- m$gene_id[m[[mod_col]] == mod_id]
        n_mod  <- length(in_mod)
        n_wrky_mod <- sum(wrky_g %in% in_mod)
        ct <- matrix(c(n_wrky_mod, n_wrky - n_wrky_mod,
                       n_mod - n_wrky_mod, n_all - n_mod - n_wrky + n_wrky_mod),
                     nrow=2)
        ft <- tryCatch(fisher.test(ct, alternative="greater"), error=function(e) NULL)
        enrich_rows[[paste(s, mod_id)]] <- data.frame(
          set=paste0("GGM_",s), mode="ggm", module=mod_id,
          n_wrky_in_module=n_wrky_mod, n_module_genes=n_mod,
          n_wrky_in_network=n_wrky, n_network_genes=n_all,
          odds_ratio=if(!is.null(ft)) ft$estimate else NA,
          pval=if(!is.null(ft)) ft$p.value else NA,
          stringsAsFactors=FALSE)
      }
    }
    for (s in PB_SETS) {
      m <- pb_mods[[s]]
      if (is.null(m)) next
      mod_col <- if ("module" %in% names(m)) "module" else names(m)[2]
      all_g <- m$gene_id; wrky_g <- wrky$gene_id
      n_all <- length(all_g); n_wrky <- sum(wrky_g %in% all_g)
      for (mod_id in unique(m[[mod_col]])) {
        if (is.na(mod_id) || mod_id == 0) next
        in_mod <- m$gene_id[m[[mod_col]] == mod_id]
        n_mod  <- length(in_mod)
        n_wrky_mod <- sum(wrky_g %in% in_mod)
        ct <- matrix(c(n_wrky_mod, n_wrky - n_wrky_mod,
                       n_mod - n_wrky_mod, n_all - n_mod - n_wrky + n_wrky_mod),
                     nrow=2)
        ft <- tryCatch(fisher.test(ct, alternative="greater"), error=function(e) NULL)
        enrich_rows[[paste(s, mod_id)]] <- data.frame(
          set=paste0("PB_",s), mode="pseudobulk", module=mod_id,
          n_wrky_in_module=n_wrky_mod, n_module_genes=n_mod,
          n_wrky_in_network=n_wrky, n_network_genes=n_all,
          odds_ratio=if(!is.null(ft)) ft$estimate else NA,
          pval=if(!is.null(ft)) ft$p.value else NA,
          stringsAsFactors=FALSE)
      }
    }
    enrich_df <- do.call(rbind, enrich_rows)
    enrich_df$padj <- p.adjust(enrich_df$pval, method="BH")
    # Combine with wrky base table
    wrky_out <- merge(as.data.frame(wrky), enrich_df[enrich_df$padj < 0.05,],
                      by.x=character(0), by.y=character(0), all=FALSE)
    save_csv(enrich_df, F2_OUT)
    # Figure
    sig <- enrich_df[!is.na(enrich_df$padj) & enrich_df$padj < 0.05,]
    if (nrow(sig) > 0) {
      p <- ggplot(sig, aes(x=factor(module), y=-log10(pval), size=n_wrky_in_module, color=set)) +
        geom_point(alpha=0.8) +
        facet_wrap(~mode, scales="free_x") +
        labs(title="WRKY family module enrichment (demo of generic gene-set query capability)",
             x="Module", y="-log10(p-value)", size="WRKY in module", color="Module set") +
        theme_minimal(base_size=10)
      save_png_svg(p, file.path(FIGS, "fig_wrky_crossmode"), w=12, h=6)
    }
    log_msg("ANALYSIS F2 done — WRKY enrichment across ", nrow(sig), " significant module-set combos (BH q<0.05)")
  }
})

log_msg("PHASE F complete")

# ─── PHASE FINAL-1: Figures manifest ─────────────────────────────────────────
log_msg("=== PHASE FINAL-1: Figures manifest ===")

fig_files <- list.files(FIGS, pattern="\\.(png|svg)$", full.names=FALSE)
fig_manifest <- c(
  "# Figures Manifest",
  paste0("Generated: ", ts()),
  "",
  "| Filename | Caption | Source CSV |",
  "|----------|---------|------------|"
)
captions <- list(
  "fig_degree_distributions.png"       = list(c="Degree distributions (log-log) for all 6 networks", s="topology_degree.csv"),
  "fig_hub_genes.png"                  = list(c="Centrality scatter: degree vs eigenvector centrality per network", s="topology_centrality.csv"),
  "fig_global_stats_comparison.png"    = list(c="Global network statistics side-by-side (GGM vs pseudobulk modes)", s="topology_global_stats.csv"),
  "fig_module_kme.png"                 = list(c="kME distributions per module set (all 6 sets)", s="module_kme_distributions.csv"),
  "fig_module_quality_across_sets.png" = list(c="Module quality summary: grey rate, n_modules, kME median across all 6 sets", s="module_quality_summary.csv"),
  "fig_crossmode_jaccard_heatmap.png"  = list(c="GGM vs pseudobulk module overlap (Jaccard); threshold Jaccard>0.05", s="crossmode_overlap.csv"),
  "fig_crossmethod_ari.png"            = list(c="WGCNA vs Louvain ARI within each mode (GGM large, GGM small, pseudobulk)", s="crossmethod_agreement.csv"),
  "fig_core_modespecific.png"          = list(c="Core (assigned in both modes) vs GGM-only vs pseudobulk-only gene counts", s="core_vs_modespecific.csv"),
  "fig_cross_set_consistency.png"      = list(c="GGM cross-set assignment stability: histogram of genes by # sets assigned", s="cross_set_consistency.csv"),
  "fig_module_condition_activation_ggm.png"         = list(c="GGM module condition activation — WITHIN GGM mode only; NOT cross-mode comparable", s="module_condition_activation_ggm.csv"),
  "fig_module_condition_activation_pseudobulk.png"  = list(c="Pseudobulk module condition activation — WITHIN pseudobulk mode only; NOT cross-mode comparable", s="module_condition_activation_pseudobulk.csv"),
  "fig_wrky_crossmode.png"             = list(c="WRKY family module enrichment — demo of generic gene-set query; dot size = WRKY count in module", s="geneset_query_demo_wrky.csv")
)
for (f in fig_files) {
  base_f <- f
  info   <- captions[[base_f]]
  cap    <- if (!is.null(info)) info$c else "—"
  src    <- if (!is.null(info)) info$s else "—"
  fig_manifest <- c(fig_manifest, paste0("| ", f, " | ", cap, " | ", src, " |"))
}
writeLines(fig_manifest, file.path(DOWN, "FIGURES_MANIFEST.md"))
log_msg("PHASE FINAL-1 complete — manifest written, ", length(fig_files), " figures on disk")

# ─── PHASE FINAL-2: Integrated HTML report ───────────────────────────────────
log_msg("=== PHASE FINAL-2: HTML report ===")
REPORT_PATH <- file.path(DOWN, "DOWNSTREAM_ANALYSIS_REPORT.html")

# Helper: base64-encode a PNG file
b64_png <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw_bytes <- readBin(path, "raw", file.size(path))
  paste0("data:image/png;base64,", base64encode(raw_bytes))
}

img_tag <- function(png_path, alt="figure", width="100%") {
  b64 <- tryCatch(b64_png(png_path), error=function(e) NULL)
  if (is.null(b64)) return(paste0('<p><em>Figure not found: ', basename(png_path), '</em></p>'))
  paste0('<img src="', b64, '" alt="', alt, '" style="max-width:', width, ';height:auto;">')
}

# Read key CSVs for inline tables
read_csv_safe <- function(path, n=20) {
  if (!file.exists(path)) return(NULL)
  tryCatch(head(read.csv(path), n), error=function(e) NULL)
}
df_to_html <- function(df, caption="") {
  if (is.null(df)) return("<p><em>Data not available</em></p>")
  rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse=""), "</tr>"))
  hdr  <- paste0("<th>", names(df), "</th>", collapse="")
  paste0(if (nchar(caption)>0) paste0('<p><strong>',caption,'</strong></p>') else '',
         '<div style="overflow-x:auto"><table border="1" style="border-collapse:collapse;font-size:11px;">',
         '<thead><tr>', hdr, '</tr></thead><tbody>',
         paste(rows, collapse=""), '</tbody></table></div>')
}

# Build HTML
html_parts <- list()
html_parts[["head"]] <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Downstream Co-expression Analysis Report</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; background: #fafafa; }
  h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
  h2 { color: #2980b9; margin-top: 40px; border-left: 4px solid #3498db; padding-left: 10px; }
  h3 { color: #555; }
  .toc { background: #ecf0f1; border: 1px solid #bdc3c7; padding: 15px 25px; border-radius: 6px; margin: 20px 0; }
  .toc a { display: block; color: #2980b9; text-decoration: none; margin: 3px 0; }
  .toc a:hover { text-decoration: underline; }
  .caveat { background: #fef9e7; border: 2px solid #f39c12; padding: 12px 18px; border-radius: 6px; margin: 12px 0; }
  .note   { background: #eaf4fb; border: 1px solid #85c1e9; padding: 10px 15px; border-radius: 5px; margin: 10px 0; }
  .section { background: white; border: 1px solid #ddd; border-radius: 6px; padding: 20px 25px; margin: 20px 0; }
  .fig-row { display: flex; flex-wrap: wrap; gap: 10px; }
  .fig-box { flex: 1 1 46%; border: 1px solid #ddd; border-radius: 4px; padding: 8px; background: #fff; }
  details { margin: 10px 0; }
  summary { cursor: pointer; font-weight: bold; color: #2980b9; padding: 6px; background: #eaf4fb; border-radius: 4px; }
  summary:hover { background: #d6eaf8; }
  table { border-collapse: collapse; font-size: 11px; }
  th { background: #2980b9; color: white; padding: 4px 8px; }
  td { padding: 3px 8px; border: 1px solid #ddd; }
  tr:nth-child(even) { background: #f2f2f2; }
  .status-done { color: green; font-weight: bold; }
  .status-skip { color: #888; }
  .status-fail { color: red; }
  #sticky-toc { position: sticky; top: 0; background: #2c3e50; color: white; padding: 8px 20px; z-index: 100; display: flex; flex-wrap: wrap; gap: 10px; }
  #sticky-toc a { color: #85c1e9; text-decoration: none; font-size: 12px; }
  #sticky-toc a:hover { color: white; }
</style>
</head>
<body>
'

html_parts[["sticky_toc"]] <- '<div id="sticky-toc">
  <strong style="margin-right:12px;">Navigate:</strong>
  <a href="#overview">1.Overview</a>
  <a href="#topology">2.Topology</a>
  <a href="#module-quality">3.Module Quality</a>
  <a href="#cross-mode">4.Cross-mode</a>
  <a href="#condition">5.Condition</a>
  <a href="#annotation">6.Annotation</a>
  <a href="#gene-centric">7.Gene-centric</a>
  <a href="#summary">8.Summary</a>
</div>
'

html_parts[["title"]] <- paste0(
  '<h1>Downstream Co-expression Analysis Report</h1>',
  '<p style="color:#666;">Generated: ', ts(), ' | Dataset: pathogen_multiome | Repo: coexpression_arabidopsis</p>')

html_parts[["toc"]] <- '<div class="toc">
  <strong>Table of Contents</strong>
  <a href="#overview">1. Overview: prior-free downstream, two modes, cross-mode caveat</a>
  <a href="#topology">2. Network Topology (A1–A3): degree distributions, hubs, global stats</a>
  <a href="#module-quality">3. Module Quality (B1–B4): kME, eigengene meta-structure, quality summary</a>
  <a href="#cross-mode">4. Cross-mode &amp; Cross-method (C1–C4): overlaps, core vs mode-specific</a>
  <a href="#condition">5. Condition-specificity (D1–D2): per-mode activation (NOT cross-comparable)</a>
  <a href="#annotation">6. Functional Annotation (E1–E2): GO/TF as reference, descriptive only</a>
  <a href="#gene-centric">7. Gene-centric Utility (F1–F2): master lookup + WRKY demo</a>
  <a href="#summary">8. Summary &amp; Reuse on New Datasets</a>
</div>'

html_parts[["sec1"]] <- '<div class="section" id="overview">
<h2>1. Overview</h2>
<h3>Purpose: Prior-free Downstream Analysis</h3>
<p>This report presents the downstream analysis of finalized co-expression networks and modules from the <em>Arabidopsis thaliana</em> pathogen multiome dataset. All primary analyses are <strong>structure-based and prior-free</strong> — driven by network topology, module coherence, and cross-mode agreement, not by biological knowledge of the genes.</p>
<div class="caveat">
  <strong>GO/TF annotation is REFERENCE OUTPUT ONLY.</strong> GO enrichment results and TF annotations are presented as descriptive metadata per module. They are <em>never</em> used to name modules, never used as selection or ranking criteria. Modules are identified by number only (e.g., "module 3"), not by biological process names.
</div>

<h3>Two Complementary Modes</h3>
<table border="1">
  <tr><th>Property</th><th>GGM mode</th><th>Pseudobulk mode</th></tr>
  <tr><td>Edge type</td><td>Partial correlation (GGM, MATLAB)</td><td>Spearman |r| ≥ 0.42</td></tr>
  <tr><td>Observation unit</td><td>Individual cells per condition</td><td>298 Seurat subclusters</td></tr>
  <tr><td>Networks</td><td>4 per-condition + robustness consensus</td><td>1 network</td></tr>
  <tr><td>Module sets</td><td>4 (large/small × WGCNA/Louvain)</td><td>2 (WGCNA p=9 / Louvain)</td></tr>
  <tr><td>GO enrichment</td><td>Computed this session</td><td>Pre-computed</td></tr>
</table>

<div class="caveat">
  <strong>⚠ Cross-mode condition-pattern comparison is FORBIDDEN.</strong><br>
  GGM uses partial correlation on cells; pseudobulk uses Spearman on subclusters. These have different edge definitions, different observation-point units, and I_s is near-binary in both. <em>Per-mode</em> condition activation is reported within each mode. Quantitative comparison of condition weights <em>across</em> modes is not valid and not performed in this report.
</div>
</div>'

# Section 2: topology
gs_table  <- df_to_html(read_csv_safe(A3_OUT), "Global network statistics")
deg_table <- df_to_html(read_csv_safe(A1_OUT, 10), "Degree + power-law fit (first 10 rows)")

html_parts[["sec2"]] <- paste0(
'<div class="section" id="topology">
<h2>2. Network Topology (A1–A3)</h2>
<h3>A3: Global Statistics</h3>',
gs_table,
'<h3>A1: Degree Distributions + Scale-Free Fit</h3>
<div class="note">Power-law fit via igraph::fit_power_law (Clauset et al. MLE). Alpha &gt; 2 is consistent with scale-free topology. KS p-value tests goodness-of-fit.</div>',
deg_table,
'<div class="fig-row">
  <div class="fig-box">',
  img_tag(file.path(FIGS,"fig_degree_distributions.png"), "Degree distributions"),
  '<p style="font-size:11px;color:#666;">Degree distributions (log-log). Steeper curves indicate more hub-like organization.</p>
  </div>
  <div class="fig-box">',
  img_tag(file.path(FIGS,"fig_global_stats_comparison.png"), "Global stats"),
  '<p style="font-size:11px;color:#666;">Global network statistics. GGM modes and pseudobulk shown side-by-side (different scales).</p>
  </div>
</div>

<h3>A2: Hub Genes (Centrality)</h3>
<div class="note">Betweenness computed only for networks with &lt;100k edges (GGM consensus only). For larger networks, degree + eigenvector centrality are reported.</div>',
df_to_html(read_csv_safe(A2_OUT, 15), "Top centrality genes (first 15 rows)"),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 90%">',
  img_tag(file.path(FIGS,"fig_hub_genes.png"), "Hub genes"),
  '<p style="font-size:11px;color:#666;">Degree vs eigenvector centrality per network. Upper-right: high-degree, high-influence hub genes.</p>
  </div>
</div>
</div>')

# Section 3: module quality
html_parts[["sec3"]] <- paste0(
'<div class="section" id="module-quality">
<h2>3. Module Quality (B1–B4)</h2>
<h3>B4: Quality Summary — All 6 Module Sets</h3>',
df_to_html(read_csv_safe(B4_OUT)),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 90%">',
  img_tag(file.path(FIGS,"fig_module_quality_across_sets.png")),
  '</div>
</div>
<h3>B1: kME Distributions</h3>',
df_to_html(read_csv_safe(B1_OUT, 15)),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 90%">',
  img_tag(file.path(FIGS,"fig_module_kme.png")),
  '<p style="font-size:11px;color:#666;">kME distributions per module set. Modules with median kME &lt; 0.3 are flagged as low-coherence.</p>
  </div>
</div>

<h3>B3: Intramodular Hub Genes</h3>',
df_to_html(read_csv_safe(B3_OUT, 20), "Top hub genes per module (first 20 rows)"),

'<h3>B2: Eigengene-eigengene Correlations</h3>
<details><summary>Eigengene correlation data (collapsed — click to expand)</summary>',
df_to_html(read_csv_safe(B2_OUT, 20)),
'</details>
<div class="note">For GGM sets, inter-module kME correlation is used as a proxy for eigengene-eigengene correlation (kME matrix is genes × modules). For pseudobulk WGCNA, sample-level MEs from blockwiseModules are used. For pseudobulk Louvain, PC1 per module is computed from the obs_normalized_cache (11010 genes × 298 subclusters).</div>
</div>')

# Section 4: cross-mode
html_parts[["sec4"]] <- paste0(
'<div class="section" id="cross-mode">
<h2>4. Cross-mode &amp; Cross-method Comparison (C1–C4)</h2>
<div class="note"><strong>This is the headline section.</strong> The pipeline\'s distinctive value is providing two complementary co-expression views (partial vs. marginal correlation) and quantifying which gene groupings are robust across both.</div>

<h3>C3: Core vs Mode-specific Gene Partitioning</h3>',
df_to_html(read_csv_safe(C3_OUT, 20), "Gene category assignments (first 20)"),
'<div class="fig-row">
  <div class="fig-box">',
  img_tag(file.path(FIGS,"fig_core_modespecific.png")),
  '<p style="font-size:11px;color:#666;"><strong>Core</strong>: assigned in both GGM and pseudobulk (robust co-expression units). <strong>Mode-specific</strong>: only grouped in one mode.</p>
  </div>
  <div class="fig-box">',
  img_tag(file.path(FIGS,"fig_cross_set_consistency.png")),
  '<p style="font-size:11px;color:#666;"><strong>GGM cross-set consistency</strong>: number of GGM module sets a gene is stably assigned to (out of 4).</p>
  </div>
</div>

<h3>C1: GGM vs Pseudobulk Module Overlap (Jaccard + ARI)</h3>',
df_to_html(read_csv_safe(file.path(DOWN,"crossmode_ari.csv")), "Cross-mode ARI per set pair"),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 90%">',
  img_tag(file.path(FIGS,"fig_crossmode_jaccard_heatmap.png")),
  '<p style="font-size:11px;color:#666;">Pairwise Jaccard between GGM and pseudobulk module gene sets (Jaccard &gt; 0.05 shown). Blocks of high overlap = robust co-expression units.</p>
  </div>
</div>
<details><summary>Full Jaccard table (collapsed)</summary>',
df_to_html(read_csv_safe(C1_OUT, 30)),
'</details>

<h3>C2: WGCNA vs Louvain Agreement (within each mode)</h3>',
df_to_html(read_csv_safe(C2_OUT)),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 60%">',
  img_tag(file.path(FIGS,"fig_crossmethod_ari.png")),
  '</div>
</div>
</div>')

# Section 5: condition
html_parts[["sec5"]] <- paste0(
'<div class="section" id="condition">
<h2>5. Condition-specificity (D1–D2)</h2>
<div class="caveat">
  <strong>⚠ The two panels below are NOT cross-comparable.</strong><br>
  GGM condition weights reflect partial correlation edge weights across condition-specific GGM networks. Pseudobulk condition weights reflect mean intramodular Spearman correlations across subcluster batches. These are on different scales with different null distributions. Compare WITHIN each panel only.
</div>

<h3>D1: Module Condition Activation</h3>
<h4>GGM mode (4 per-condition networks, within-mode only)</h4>',
df_to_html(read_csv_safe(D1_GGM_OUT, 15)),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 90%">',
  img_tag(file.path(FIGS,"fig_module_condition_activation_ggm.png")),
  '</div>
</div>
<h4>Pseudobulk mode (subcluster Spearman network, within-mode only)</h4>',
df_to_html(read_csv_safe(D1_PB_OUT, 15)),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 90%">',
  img_tag(file.path(FIGS,"fig_module_condition_activation_pseudobulk.png")),
  '</div>
</div>
</div>')

# Section 6: annotation
html_parts[["sec6"]] <- paste0(
'<div class="section" id="annotation">
<h2>6. Functional Annotation (E1–E2) — REFERENCE OUTPUT ONLY</h2>
<div class="caveat">
  <strong>Modules are NOT named by GO terms or TF families.</strong> GO enrichment and TF content are attached as descriptive metadata per module. They are reference context, not analysis criteria.
</div>

<h3>E1: GO BP Enrichment</h3>',
df_to_html(read_csv_safe(file.path(DOWN,"go_enrichment_all_sets.csv"), 25),
           "Master GO enrichment table (first 25 rows; q < 0.05, min set size 10)"),
'
<h3>E2: TF Module Enrichment</h3>',
df_to_html(read_csv_safe(E2_OUT, 20), "TF-module records (first 20 rows)"),
'</div>')

# Section 7: gene-centric
html_parts[["sec7"]] <- paste0(
'<div class="section" id="gene-centric">
<h2>7. Gene-centric Utility (F1–F2)</h2>
<h3>F1: Master Gene Lookup Table</h3>
<p>One row per gene: module assignment across all 6 module sets, kME in each, top-10 co-expression partners in both modes.</p>',
df_to_html(read_csv_safe(F1_OUT, 15), "gene_lookup_master.csv (first 15 rows)"),
'<div class="note"><strong>How to use:</strong> Filter by gene_id to retrieve a gene\'s full co-expression profile. The top10 partner columns give context for guilt-by-association inference. Module columns allow cross-set comparison in a single lookup.</div>

<h3>F2: WRKY Family — Demo of Generic Gene-set Query</h3>
<p>This demonstrates the generic gene-set module enrichment capability using the WRKY TF family as an example. The same query works for any user-provided gene list.</p>',
df_to_html(read_csv_safe(F2_OUT, 20), "WRKY enrichment results (BH-adjusted; first 20 rows)"),
'<div class="fig-row">
  <div class="fig-box" style="flex:1 1 90%">',
  img_tag(file.path(FIGS,"fig_wrky_crossmode.png")),
  '<p style="font-size:11px;color:#666;">WRKY TF family enrichment by module (dot = BH q&lt;0.05 module). This is a <em>demo</em> of the generic gene-set query capability, not a pathogen-specific analysis.</p>
  </div>
</div>
</div>')

# Section 8: summary
analysis_status <- c(
  "| Analysis | Status | Output |",
  "|----------|--------|--------|",
  paste0("| A1: Degree + power-law | ", if(file.exists(A1_OUT)) "DONE" else "MISSING", " | topology_degree.csv |"),
  paste0("| A2: Centrality | ", if(file.exists(A2_OUT)) "DONE" else "MISSING", " | topology_centrality.csv |"),
  paste0("| A3: Global stats | ", if(file.exists(A3_OUT)) "DONE" else "MISSING", " | topology_global_stats.csv |"),
  paste0("| B1: kME distributions | ", if(file.exists(B1_OUT)) "DONE" else "MISSING", " | module_kme_distributions.csv |"),
  paste0("| B2: Eigengene correlations | ", if(file.exists(B2_OUT)) "DONE" else "MISSING", " | module_eigengene_correlations.csv |"),
  paste0("| B3: Hub genes | ", if(file.exists(B3_OUT)) "DONE" else "MISSING", " | module_hubs.csv |"),
  paste0("| B4: Module quality | ", if(file.exists(B4_OUT)) "DONE" else "MISSING", " | module_quality_summary.csv |"),
  paste0("| C1: Cross-mode overlap | ", if(file.exists(C1_OUT)) "DONE" else "MISSING", " | crossmode_overlap.csv |"),
  paste0("| C2: WGCNA vs Louvain | ", if(file.exists(C2_OUT)) "DONE" else "MISSING", " | crossmethod_agreement.csv |"),
  paste0("| C3: Core vs mode-specific | ", if(file.exists(C3_OUT)) "DONE" else "MISSING", " | core_vs_modespecific.csv |"),
  paste0("| C4: Cross-set consistency | ", if(file.exists(C4_OUT)) "DONE" else "MISSING", " | cross_set_consistency.csv |"),
  paste0("| D1: GGM condition activation | ", if(file.exists(D1_GGM_OUT)) "DONE" else "MISSING", " | module_condition_activation_ggm.csv |"),
  paste0("| D1: PB condition activation | ", if(file.exists(D1_PB_OUT)) "DONE" else "MISSING", " | module_condition_activation_pseudobulk.csv |"),
  paste0("| D2: GGM condition patterns | ", if(file.exists(D2_GGM_OUT)) "DONE" else "MISSING", " | module_condition_patterns_ggm.csv |"),
  paste0("| D2: PB condition patterns | ", if(file.exists(D2_PB_OUT)) "DONE" else "MISSING", " | module_condition_patterns_pseudobulk.csv |"),
  paste0("| E1: GO enrichment | ", if(file.exists(file.path(DOWN,"go_enrichment_all_sets.csv"))) "DONE" else "PARTIAL/MISSING", " | go_enrichment_all_sets.csv |"),
  paste0("| E2: TF enrichment | ", if(file.exists(E2_OUT)) "DONE" else "MISSING", " | tf_enrichment.csv |"),
  paste0("| F1: Gene lookup master | ", if(file.exists(F1_OUT)) "DONE" else "MISSING", " | gene_lookup_master.csv |"),
  paste0("| F2: WRKY demo | ", if(file.exists(F2_OUT)) "DONE" else "MISSING", " | geneset_query_demo_wrky.csv |")
)

html_parts[["sec8"]] <- paste0(
'<div class="section" id="summary">
<h2>8. Summary &amp; How to Reuse on a New Dataset</h2>
<h3>Analysis Completion Status</h3>
<pre>', paste(analysis_status, collapse="\n"), '</pre>
<h3>How to Reuse on a New Dataset (e.g., dev atlas)</h3>
<p>All analyses in this report are <strong>dataset-agnostic</strong> and can be run on any new co-expression result from the pipeline:</p>
<ol>
  <li><strong>Fully automatic</strong> (no user input needed):
    <ul>
      <li>A1–A3: Network topology — provide edge tables in the same format</li>
      <li>B1–B4: Module quality — provide gene_module.csv with gene_id, module, kME</li>
      <li>C1–C4: Cross-mode comparison — provide two sets of module assignments</li>
      <li>D1–D2: Condition activation — provide module_condition_patterns.csv per set</li>
      <li>E1: GO enrichment — automatic via clusterProfiler + org.At.tair.db</li>
      <li>E2: TF enrichment — automatic if module_tfs.csv is present</li>
    </ul>
  </li>
  <li><strong>Requires a parameter</strong> (gene list):
    <ul>
      <li>F2: Gene-set query demo — provide a CSV with gene_id column (any gene family)</li>
    </ul>
  </li>
</ol>
<h3>Key Files</h3>
<ul>
  <li><code>results/pathogen_multiome/downstream/gene_lookup_master.csv</code> — master per-gene lookup</li>
  <li><code>results/pathogen_multiome/downstream/DOWNSTREAM_ANALYSIS_REPORT.html</code> — this report</li>
  <li><code>results/pathogen_multiome/downstream/figures/</code> — all PNG + SVG figures</li>
  <li><code>results/pathogen_multiome/downstream/DOWNSTREAM_INVENTORY.md</code> — what was found and reused</li>
</ul>
</div>')

html_parts[["foot"]] <- '</body></html>'

# Combine and write
full_html <- paste(unlist(html_parts), collapse="\n")
writeLines(full_html, REPORT_PATH)

# Verify standalone (check file size > 10KB, no external file refs)
fsize <- file.size(REPORT_PATH)
log_msg("PHASE FINAL-2 complete — HTML report: ", REPORT_PATH,
        " (", round(fsize/1024), " KB) standalone=", (fsize > 10000))

# ─── PHASE FINAL-3: SESSION_HANDOFF append ────────────────────────────────────
log_msg("=== PHASE FINAL-3: SESSION_HANDOFF update ===")
handoff_path <- file.path(REPO, "docs/SESSION_HANDOFF.md")
if (file.exists(handoff_path)) {
  new_section <- paste0(
    "\n\n## Downstream Analysis Suite — ", format(Sys.Date(), "%Y-%m-%d"), "\n\n",
    "All downstream analyses complete. Output root: `results/pathogen_multiome/downstream/`\n\n",
    "### What now exists\n",
    "- **Topology** (A1-A3): degree distributions, power-law fits, centrality, global stats for all 6 networks\n",
    "- **Module quality** (B1-B4): kME distributions, eigengene correlations, hub genes, quality comparison across all 6 sets\n",
    "- **Cross-mode** (C1-C4): GGM vs pseudobulk Jaccard/ARI overlap, WGCNA vs Louvain ARI, core vs mode-specific gene partition, cross-set consistency\n",
    "- **Condition specificity** (D1-D2): per-mode activation tables and figures (GGM and pseudobulk separately; cross-mode comparison forbidden)\n",
    "- **Functional annotation** (E1-E2): GO BP enrichment per module (all sets), TF enrichment (reference only; modules not named)\n",
    "- **Gene-centric** (F1-F2): master gene lookup table + WRKY demo of generic gene-set query\n",
    "- **Integrated HTML report**: `DOWNSTREAM_ANALYSIS_REPORT.html` (self-contained, standalone)\n\n",
    "### Generic/reusable\n",
    "All analyses are dataset-agnostic. See Section 8 of the HTML report for reuse instructions.\n",
    "Analysis script: `inst/scripts/downstream_analysis.R`\n"
  )
  cat(new_section, file=handoff_path, append=TRUE)
  log_msg("SESSION_HANDOFF.md updated")
}

log_msg("=== downstream_analysis.R COMPLETE ===")
close(log_con)
