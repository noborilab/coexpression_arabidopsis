#!/usr/bin/env Rscript
# Overnight analysis: core/mode-specific gene characterization +
# marker-vs-co-expression co-detection benchmark + PRIMER case study.
# Phases 0–7: runs end-to-end unattended; saves incrementally; one failure
# never aborts the run.
# Usage: Rscript inst/scripts/core_marker_analysis.R

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# ── Install missing packages silently ────────────────────────────────────────
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
  requireNamespace(pkg, quietly = TRUE)
}
ensure_pkg("base64enc")

# ── Logging ───────────────────────────────────────────────────────────────────
LOG_PATH <- "logs/core_marker_overnight.log"
dir.create("logs", showWarnings = FALSE)
log_con  <- file(LOG_PATH, open = "at")
logmsg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
  message(msg)
  cat(msg, "\n", file = log_con, append = TRUE)
}
logmsg("=== core_marker_analysis.R start ===")

# ── Paths ─────────────────────────────────────────────────────────────────────
DROPBOX  <- paste0("/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/",
                   "SALK_clowd/Projects/SA_PTI_ETI_single_cell/",
                   "SA_039_94_multiome_revision_rep2_9h_only")

SEURAT_PATH  <- file.path(DROPBOX, "out/_seurat_object/motifFixed/combined_filtered.rds")
MARKER_MAJOR <- file.path(DROPBOX, "out/2_clustering/RNA/markers_RNA.txt")
MARKER_IMMUNE <- file.path(DROPBOX, "out/5_subclustering/3_7_11/markers_ch_removed.txt")

RES_ROOT  <- "results/pathogen_multiome"
OUT_ROOT  <- file.path(RES_ROOT, "core_marker")
FIG_DIR   <- file.path(OUT_ROOT, "figures")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

SYMBOL_MAP_PATH <- file.path(RES_ROOT, "symbol_map.csv")
CORE_CSV        <- file.path(RES_ROOT, "downstream/core_vs_modespecific.csv")

GGM_SETS <- list(
  GGM_large_wgcna  = file.path(RES_ROOT, "official_modules/large_wgcna/gene_module.csv"),
  GGM_large_louvain = file.path(RES_ROOT, "official_modules/large_louvain/gene_module.csv"),
  GGM_small_wgcna  = file.path(RES_ROOT, "official_modules/small_wgcna/gene_module.csv"),
  GGM_small_louvain = file.path(RES_ROOT, "official_modules/small_louvain/gene_module.csv")
)
PB_SETS <- list(
  PB_wgcna   = file.path(RES_ROOT, "pseudobulk_zscore_spearman/modules_official/wgcna/module_membership.csv"),
  PB_louvain = file.path(RES_ROOT, "pseudobulk_zscore_spearman/modules_official/louvain/module_membership.csv")
)

# ── Source co-detection engine ────────────────────────────────────────────────
source("inst/scripts/codetection_eval.R")

# ── Phase status tracker ──────────────────────────────────────────────────────
phase_status <- list()
mark_phase <- function(name, status) {
  phase_status[[name]] <<- status
  logmsg(sprintf("PHASE %s: %s", name, status))
}

skip_if_exists <- function(path) {
  if (file.exists(path)) {
    logmsg(sprintf("  SKIP (exists): %s", basename(path)))
    return(TRUE)
  }
  FALSE
}

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 0: Confirm inputs + build display labels
# ─────────────────────────────────────────────────────────────────────────────
logmsg("=== PHASE 0: Confirm inputs ===")
tryCatch({

  # Load symbol map
  symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)
  sym2id  <- setNames(symbol_map$gene_id,     symbol_map$gene_symbol)
  id2sym  <- setNames(symbol_map$gene_symbol, symbol_map$gene_id)
  logmsg(sprintf("  symbol_map: %d entries", nrow(symbol_map)))

  make_label <- function(atid) {
    sym <- id2sym[atid]
    if (!is.na(sym) && nchar(trimws(sym)) > 0) paste0(sym, " (", atid, ")") else atid
  }

  # Map gene name → AT-ID (handles symbols AND AT-IDs)
  resolve_to_atid <- function(gname) {
    gname <- trimws(gname)
    if (grepl("^AT[0-9MC]G[0-9]+", gname, ignore.case = TRUE))
      return(toupper(gname))
    hit <- sym2id[gname]
    if (!is.na(hit)) return(hit)
    NA_character_
  }

  # Load Seurat object
  logmsg("  Loading Seurat object...")
  suppressPackageStartupMessages(library(Seurat))
  seurat_obj <- readRDS(SEURAT_PATH)
  n_cells <- ncol(seurat_obj)
  n_genes <- nrow(seurat_obj[["RNA"]])
  logmsg(sprintf("  Seurat: %d cells, %d genes", n_cells, n_genes))
  DefaultAssay(seurat_obj) <- "RNA"
  RNA_GENES <- rownames(seurat_obj[["RNA"]])

  # Build AT-ID → RNA rowname lookup
  atid_to_rna <- character(0)
  for (i in seq_len(nrow(symbol_map))) {
    gid  <- symbol_map$gene_id[i]
    gsym <- symbol_map$gene_symbol[i]
    if (!is.na(gsym) && gsym %in% RNA_GENES) atid_to_rna[gid] <- gsym
    else if (gid %in% RNA_GENES)             atid_to_rna[gid] <- gid
  }
  for (g in RNA_GENES[grepl("^AT[0-9MC]G[0-9]+", RNA_GENES, ignore.case = TRUE)]) {
    gup <- toupper(g)
    if (!(gup %in% names(atid_to_rna))) atid_to_rna[gup] <- g
  }
  logmsg(sprintf("  AT-ID→RNA map: %d entries (out of %d RNA genes)", length(atid_to_rna), n_genes))

  resolve_grp <- function(atids) {
    rna_names <- atid_to_rna[atids]
    rna_names[!is.na(rna_names)]
  }

  # Detect metadata columns (major cluster + subcluster)
  md_cols <- colnames(seurat_obj@meta.data)
  logmsg(sprintf("  Metadata columns: %s", paste(md_cols, collapse=", ")))

  # Read major cluster markers
  logmsg("  Reading major cluster markers...")
  mk_major_raw <- read.table(MARKER_MAJOR, sep = "\t", header = TRUE,
                              row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  logmsg(sprintf("  markers_RNA: %d rows, clusters: %s",
    nrow(mk_major_raw), paste(sort(unique(mk_major_raw$cluster)), collapse = ",")))

  # Read immune subcluster markers
  logmsg("  Reading immune subcluster markers...")
  mk_immune_raw <- read.table(MARKER_IMMUNE, sep = "\t", header = TRUE,
                               fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # Detect which column is the gene name and which is cluster
  # First column after reading with fill=TRUE may be gene name (row names aren't auto-detected)
  # Inspect column names
  logmsg(sprintf("  markers_ch_removed cols: %s", paste(colnames(mk_immune_raw), collapse=", ")))
  # Standard format: row is gene, cols include 'cluster' and 'gene'
  if ("gene" %in% colnames(mk_immune_raw) && "cluster" %in% colnames(mk_immune_raw)) {
    mk_immune_raw <- mk_immune_raw[!is.na(mk_immune_raw$cluster) &
                                     mk_immune_raw$cluster != "cluster", ]
    mk_immune_raw$cluster <- suppressWarnings(as.integer(mk_immune_raw$cluster))
    mk_immune_raw <- mk_immune_raw[!is.na(mk_immune_raw$cluster), ]
  } else {
    # Fallback: try reading with row.names=1
    mk_immune_raw <- read.table(MARKER_IMMUNE, sep = "\t", header = TRUE,
                                 row.names = 1, fill = TRUE, check.names = FALSE,
                                 stringsAsFactors = FALSE)
    mk_immune_raw$gene <- rownames(mk_immune_raw)
    mk_immune_raw$cluster <- suppressWarnings(as.integer(mk_immune_raw$cluster))
    mk_immune_raw <- mk_immune_raw[!is.na(mk_immune_raw$cluster), ]
  }
  logmsg(sprintf("  markers_ch_removed: %d rows, clusters: %s",
    nrow(mk_immune_raw), paste(sort(unique(mk_immune_raw$cluster)), collapse = ",")))

  # Build marker group lists: cluster → AT-IDs
  build_marker_groups <- function(df, src_label) {
    clusters <- sort(unique(df$cluster))
    grps <- lapply(clusters, function(cl) {
      genes_raw <- df$gene[df$cluster == cl]
      atids <- sapply(genes_raw, resolve_to_atid, USE.NAMES = FALSE)
      atids <- unique(atids[!is.na(atids)])
      atids
    })
    names(grps) <- paste0(src_label, "_c", clusters)
    grps
  }

  marker_groups_major  <- build_marker_groups(mk_major_raw,  "MajorClust")
  marker_groups_immune <- build_marker_groups(mk_immune_raw, "ImmuneSubclust")

  logmsg(sprintf("  Major marker groups: %d (sizes: %s)",
    length(marker_groups_major),
    paste(range(sapply(marker_groups_major, length)), collapse = "–")))
  logmsg(sprintf("  Immune marker groups: %d (sizes: %s)",
    length(marker_groups_immune),
    paste(range(sapply(marker_groups_immune, length)), collapse = "–")))

  # Marker multiplicity: how many clusters mark each gene?
  all_major_genes <- unlist(marker_groups_major, use.names = FALSE)
  mult_tbl <- table(all_major_genes)
  mult_df  <- data.frame(
    gene_id     = names(mult_tbl),
    n_clusters  = as.integer(mult_tbl),
    display_label = sapply(names(mult_tbl), make_label),
    stringsAsFactors = FALSE
  )
  write.csv(mult_df, file.path(OUT_ROOT, "marker_gene_multiplicity.csv"),
            row.names = FALSE, quote = TRUE)
  logmsg(sprintf("  Marker multiplicity: %d unique major-marker genes; multi-cluster: %d (%.1f%%)",
    nrow(mult_df), sum(mult_df$n_clusters > 1),
    100 * mean(mult_df$n_clusters > 1)))

  # Figure: multiplicity histogram
  p_mult <- ggplot(mult_df, aes(x = factor(pmin(n_clusters, 5)))) +
    geom_bar(fill = "#4575b4") +
    scale_x_discrete(labels = c("1","2","3","4","5+")) +
    labs(x = "Number of clusters gene marks", y = "Gene count",
         title = "Major-cluster marker multiplicity") +
    theme_minimal(base_size = 12)
  ggsave(file.path(FIG_DIR, "fig_marker_multiplicity.png"), p_mult,
         width = 5, height = 4, dpi = 300, bg = "white")

  # Load co-expression module gene lists
  load_ggm_module <- function(path, src_name) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    mods <- split(df$gene_id, df$top_module)
    names(mods) <- paste0(src_name, "_m", names(mods))
    mods
  }
  load_pb_module <- function(path, src_name) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    # Exclude grey/unassigned (module 0 for wgcna; check for 0)
    df <- df[df$module != 0, ]
    mods <- split(df$gene_id, df$module)
    names(mods) <- paste0(src_name, "_m", names(mods))
    mods
  }

  coexpr_modules <- list()
  for (nm in names(GGM_SETS)) {
    coexpr_modules <- c(coexpr_modules, load_ggm_module(GGM_SETS[[nm]], nm))
  }
  for (nm in names(PB_SETS)) {
    coexpr_modules <- c(coexpr_modules, load_pb_module(PB_SETS[[nm]], nm))
  }
  logmsg(sprintf("  Co-expression module groups: %d total", length(coexpr_modules)))

  # Sparse matrix extraction
  logmsg("  Extracting sparse count and lognorm matrices...")
  counts_mat  <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  lognorm_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  logmsg(sprintf("  counts_mat: %d × %d (class: %s)",
    nrow(counts_mat), ncol(counts_mat), class(counts_mat)[1]))

  # Detection rates for all genes
  logmsg("  Computing detection rates...")
  det_rate_all <- compute_det_rates(counts_mat)
  logmsg(sprintf("  det_rate_all: %d genes, mean=%.4f", length(det_rate_all), mean(det_rate_all)))

  # Load core/mode-specific
  core_df <- read.csv(CORE_CSV, stringsAsFactors = FALSE)
  core_genes    <- core_df$gene_id[core_df$category == "core"]
  ggm_only      <- core_df$gene_id[core_df$category == "GGM_specific"]
  pb_only        <- core_df$gene_id[core_df$category == "PB_specific"]
  logmsg(sprintf("  Core: %d, GGM-only: %d, PB-only: %d",
    length(core_genes), length(ggm_only), length(pb_only)))

  mark_phase("0", "DONE")
}, error = function(e) {
  logmsg(sprintf("PHASE 0 ERROR: %s", conditionMessage(e)))
  mark_phase("0", "FAILED")
  stop(e)
})

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 2: Benchmark — marker groups vs co-expression modules
# ─────────────────────────────────────────────────────────────────────────────
logmsg("=== PHASE 2: Co-detection benchmark ===")
bench_csv <- file.path(OUT_ROOT, "marker_vs_coexpr_codetection.csv")

tryCatch({
  if (!skip_if_exists(bench_csv)) {
    # Resolve all groups to RNA rownames
    resolve_all <- function(grp_list) {
      lapply(grp_list, resolve_grp)
    }

    rna_major  <- resolve_all(marker_groups_major)
    rna_immune <- resolve_all(marker_groups_immune)
    rna_coexpr <- resolve_all(coexpr_modules)

    # Score marker groups
    logmsg("  Scoring major cluster marker groups...")
    sc_major <- score_group_list(rna_major, counts_mat, lognorm_mat, det_rate_all,
                                  n_pair_cap = 5000L, seed = 98L, timeout_secs = 600L,
                                  log_fn = logmsg)
    df_major <- summarise_scores(sc_major, source_val = "marker_major")

    logmsg("  Scoring immune subcluster marker groups...")
    sc_immune <- score_group_list(rna_immune, counts_mat, lognorm_mat, det_rate_all,
                                   n_pair_cap = 5000L, seed = 98L, timeout_secs = 600L,
                                   log_fn = logmsg)
    df_immune <- summarise_scores(sc_immune, source_val = "marker_immune")

    logmsg("  Scoring co-expression module groups...")
    sc_coexpr <- score_group_list(rna_coexpr, counts_mat, lognorm_mat, det_rate_all,
                                   n_pair_cap = 5000L, seed = 98L, timeout_secs = 600L,
                                   log_fn = logmsg)
    df_coexpr <- summarise_scores(sc_coexpr, source_val = "coexpression")

    bench_df <- rbind(df_major, df_immune, df_coexpr)
    write.csv(bench_df, bench_csv, row.names = FALSE, quote = TRUE)
    logmsg(sprintf("  Saved benchmark: %d rows → %s", nrow(bench_df), bench_csv))
  } else {
    bench_df <- read.csv(bench_csv, stringsAsFactors = FALSE)
    # Reconstruct per-source score lists from saved CSV (limited, but adequate for summaries)
    sc_major  <- NULL; sc_immune <- NULL; sc_coexpr <- NULL
  }

  # Summary analysis
  logmsg("  Summarizing benchmark comparison...")

  bench_df$source_label <- dplyr::recode(bench_df$source,
    "marker_major"   = "Major-cluster markers",
    "marker_immune"  = "Immune-subcluster markers",
    "coexpression"   = "Co-expression modules")

  # Distribution of gap_codet by source
  bench_summary <- bench_df %>%
    group_by(source_label) %>%
    summarise(
      n_groups        = n(),
      mean_gap_codet  = mean(gap_codet, na.rm = TRUE),
      median_gap_codet = median(gap_codet, na.rm = TRUE),
      sd_gap_codet    = sd(gap_codet, na.rm = TRUE),
      mean_es_codet   = mean(effect_size_codet, na.rm = TRUE),
      mean_gap_spear  = mean(gap_spear, na.rm = TRUE),
      mean_es_spear   = mean(effect_size_spear, na.rm = TRUE),
      mean_n_genes    = mean(n_genes),
      .groups = "drop"
    )
  logmsg("  Benchmark summary:")
  for (i in seq_len(nrow(bench_summary))) {
    logmsg(sprintf("    %s: n=%d, codet_gap=%.4f (ES=%.2f), spear_gap=%.4f",
      bench_summary$source_label[i], bench_summary$n_groups[i],
      bench_summary$mean_gap_codet[i], bench_summary$mean_es_codet[i],
      bench_summary$mean_gap_spear[i]))
  }
  write.csv(bench_summary, file.path(OUT_ROOT, "benchmark_summary.csv"),
            row.names = FALSE, quote = TRUE)

  # Marker-weakness #1: multi-cluster markers
  mult_df <- read.csv(file.path(OUT_ROOT, "marker_gene_multiplicity.csv"), stringsAsFactors = FALSE)
  multi_genes  <- mult_df$gene_id[mult_df$n_clusters > 1]
  single_genes <- mult_df$gene_id[mult_df$n_clusters == 1]

  # Annotate major marker pairs
  if (!is.null(sc_major)) {
    multi_gap <- lapply(sc_major, function(r) {
      if (is.null(r)) return(NULL)
      g1_multi <- r$g1_ids %in% multi_genes
      g2_multi <- r$g2_ids %in% multi_genes
      has_multi <- g1_multi | g2_multi
      data.frame(
        group_id    = r$group_id,
        pair_type   = ifelse(has_multi, "has_multi_marker", "single_marker_only"),
        codet_gap   = r$codet_pairs - r$null_codet_pairs,
        stringsAsFactors = FALSE
      )
    })
    multi_gap_df <- do.call(rbind, multi_gap[!sapply(multi_gap, is.null)])
    write.csv(multi_gap_df, file.path(OUT_ROOT, "multicluster_marker_pairs.csv"),
              row.names = FALSE, quote = TRUE)
    mg_sum <- multi_gap_df %>%
      group_by(pair_type) %>%
      summarise(mean_codet_gap = mean(codet_gap, na.rm=TRUE),
                n_pairs = n(), .groups = "drop")
    logmsg("  Multi-cluster marker penalty:")
    for (i in seq_len(nrow(mg_sum))) {
      logmsg(sprintf("    %s: mean_gap=%.4f (n=%d pairs)", mg_sum$pair_type[i],
        mg_sum$mean_codet_gap[i], mg_sum$n_pairs[i]))
    }
  }

  # Figure: co-detection gap distribution by source
  bench_df_plot <- bench_df[!is.na(bench_df$gap_codet), ]
  p_bench <- ggplot(bench_df_plot, aes(x = source_label, y = gap_codet, fill = source_label)) +
    geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    labs(x = NULL, y = "Co-detection gap (within − matched null)",
         title = "Single-cell co-detection: marker groups vs co-expression modules") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(FIG_DIR, "fig_codetection_by_source.png"), p_bench,
         width = 7, height = 5, dpi = 300, bg = "white")

  # Figure: granularity — major vs immune subcluster markers
  bench_marker <- bench_df[bench_df$source %in% c("marker_major", "marker_immune"), ]
  p_gran <- ggplot(bench_marker, aes(x = source_label, y = gap_codet, fill = source_label)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    labs(x = NULL, y = "Co-detection gap",
         title = "Marker granularity effect: major vs immune-subcluster") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(FIG_DIR, "fig_granularity_effect.png"), p_gran,
         width = 5, height = 4, dpi = 300, bg = "white")

  # Write BENCHMARK_FINDINGS.md
  bmf_lines <- c(
    "# BENCHMARK FINDINGS",
    "",
    "## Central question",
    "Do co-expression modules show stronger single-cell co-detection than",
    "marker-gene groups, on the same matrix, vs expression-frequency-matched null?",
    "",
    "## Results",
    "",
    sprintf("Total groups evaluated: %d (%d major markers, %d immune markers, %d co-expression modules)",
      nrow(bench_df),
      sum(bench_df$source == "marker_major", na.rm=TRUE),
      sum(bench_df$source == "marker_immune", na.rm=TRUE),
      sum(bench_df$source == "coexpression", na.rm=TRUE)),
    "",
    "### Co-detection gap by source (within − matched-null, higher = tighter co-detection)"
  )
  for (i in seq_len(nrow(bench_summary))) {
    bmf_lines <- c(bmf_lines, sprintf(
      "- %s: mean gap = %.4f (ES = %.2f), %d groups",
      bench_summary$source_label[i], bench_summary$mean_gap_codet[i],
      bench_summary$mean_es_codet[i], bench_summary$n_groups[i]))
  }
  bmf_lines <- c(bmf_lines, "",
    "### Marker weakness #1: multi-cluster markers",
    "Genes marking >1 cluster are not exclusive co-expression units.",
    "Quantified by comparing co-detection gap for pairs involving multi-cluster",
    "markers vs single-cluster markers (see multicluster_marker_pairs.csv).",
    "",
    "### Confounds controlled",
    "- Same single-cell matrix used for all groups.",
    "- Matched-null draws pairs of genes with the SAME detection-rate decile as the",
    "  real pair. This removes the trivial confound that highly-expressed genes",
    "  co-detect at high rates by chance.",
    "- Co-expression modules differ from marker groups in size; effect sizes (not raw",
    "  gaps) are the appropriate comparison.",
    "",
    sprintf("Generated: %s", format(Sys.time()))
  )
  writeLines(bmf_lines, file.path(OUT_ROOT, "BENCHMARK_FINDINGS.md"))

  mark_phase("2", "DONE")
}, error = function(e) {
  logmsg(sprintf("PHASE 2 ERROR: %s", conditionMessage(e)))
  mark_phase("2", paste0("FAILED: ", conditionMessage(e)))
})

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3: Core and mode-specific gene analysis
# ─────────────────────────────────────────────────────────────────────────────
logmsg("=== PHASE 3: Core/mode-specific analysis ===")
core_char_csv <- file.path(OUT_ROOT, "core_modespecific_characterization.csv")

tryCatch({
  # 3a: Co-detection for core / GGM-only / PB-only
  if (!skip_if_exists(core_char_csv)) {
    set_groups <- list(
      core     = resolve_grp(core_genes),
      GGM_only = resolve_grp(ggm_only),
      PB_only  = resolve_grp(pb_only)
    )
    logmsg("  Scoring core/mode-specific groups...")
    sc_sets <- score_group_list(set_groups, counts_mat, lognorm_mat, det_rate_all,
                                 n_pair_cap = 5000L, seed = 98L, timeout_secs = 600L,
                                 log_fn = logmsg)

    # 3b: Expression character per gene per set
    all_set_genes <- c(
      setNames(rep("core",     length(core_genes)), core_genes),
      setNames(rep("GGM_only", length(ggm_only)),   ggm_only),
      setNames(rep("PB_only",  length(pb_only)),     pb_only)
    )
    char_rows <- lapply(names(all_set_genes), function(atid) {
      rna_nm <- atid_to_rna[atid]
      if (is.na(rna_nm)) return(NULL)
      dr  <- det_rate_all[rna_nm]
      mn  <- mean(as.numeric(lognorm_mat[rna_nm, ]), na.rm = TRUE)
      data.frame(gene_id = atid, display_label = make_label(atid),
                 set = all_set_genes[[atid]], det_rate = dr,
                 mean_lognorm = mn, stringsAsFactors = FALSE)
    })
    char_df <- do.call(rbind, char_rows[!sapply(char_rows, is.null)])
    write.csv(char_df, core_char_csv, row.names = FALSE, quote = TRUE)
    logmsg(sprintf("  Saved core_char: %d genes", nrow(char_df)))
  } else {
    char_df <- read.csv(core_char_csv, stringsAsFactors = FALSE)
    sc_sets <- NULL
  }

  # 3c: Marker overlap (hypergeometric)
  all_mk_atids  <- unique(unlist(marker_groups_major, use.names = FALSE))
  set_labels <- c("core", "GGM_only", "PB_only")
  set_gene_lists <- list(
    core     = core_genes,
    GGM_only = ggm_only,
    PB_only  = pb_only
  )
  N_bg <- length(unique(c(names(det_rate_all), all_mk_atids)))  # background

  marker_overlap_rows <- lapply(set_labels, function(sl) {
    set_g <- set_gene_lists[[sl]]
    K <- length(set_g)
    M <- length(all_mk_atids)
    n_overlap <- sum(set_g %in% all_mk_atids)
    # Hypergeometric: prob of >= n_overlap by chance
    pval <- phyper(n_overlap - 1, M, N_bg - M, K, lower.tail = FALSE)
    data.frame(set = sl, n_genes = K, n_marker_genes = M,
               n_overlap = n_overlap, expected = K * M / N_bg,
               enrichment_ratio = n_overlap / (K * M / N_bg),
               p_hypergeometric = pval,
               stringsAsFactors = FALSE)
  })
  overlap_df <- do.call(rbind, marker_overlap_rows)
  write.csv(overlap_df, file.path(OUT_ROOT, "core_marker_overlap.csv"),
            row.names = FALSE, quote = TRUE)
  logmsg("  Marker-set overlap:")
  for (i in seq_len(nrow(overlap_df)))
    logmsg(sprintf("    %s: overlap=%d, ratio=%.2f, p=%.2e",
      overlap_df$set[i], overlap_df$n_overlap[i],
      overlap_df$enrichment_ratio[i], overlap_df$p_hypergeometric[i]))

  # Figures
  char_df$set <- factor(char_df$set, levels = set_labels)
  p_expr <- ggplot(char_df, aes(x = set, y = det_rate, fill = set)) +
    geom_violin(alpha = 0.6) + geom_boxplot(width = 0.1, alpha = 0.8, outlier.size = 0.5) +
    labs(x = NULL, y = "Per-cell detection rate",
         title = "Expression character: core vs mode-specific genes") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(FIG_DIR, "fig_core_vs_modespecific_expression.png"), p_expr,
         width = 5, height = 4, dpi = 300, bg = "white")

  p_ovlp <- ggplot(overlap_df, aes(x = set, y = enrichment_ratio, fill = set)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = NULL, y = "Enrichment ratio (vs background)",
         title = "Marker gene overlap: core vs mode-specific") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(FIG_DIR, "fig_core_marker_overlap.png"), p_ovlp,
         width = 5, height = 4, dpi = 300, bg = "white")

  # Co-detection figure for core/mode-specific
  if (!is.null(sc_sets)) {
    sc_sum_rows <- lapply(names(sc_sets), function(nm) {
      r <- sc_sets[[nm]]
      if (is.null(r)) return(NULL)
      data.frame(set = nm,
                 codet_gap = r$codet_pairs - r$null_codet_pairs,
                 stringsAsFactors = FALSE)
    })
    sc_sum_df <- do.call(rbind, sc_sum_rows[!sapply(sc_sum_rows, is.null)])
    p_core_co <- ggplot(sc_sum_df, aes(x = set, y = codet_gap, fill = set)) +
      geom_violin(alpha = 0.6, draw_quantiles = c(0.5)) +
      geom_boxplot(width = 0.1, alpha = 0.8) +
      labs(x = NULL, y = "Co-detection gap (within − null)",
           title = "Co-detection: core vs mode-specific gene sets") +
      theme_minimal(base_size = 11) + theme(legend.position = "none")
    ggsave(file.path(FIG_DIR, "fig_core_codetection.png"), p_core_co,
           width = 5, height = 4, dpi = 300, bg = "white")
  }

  mark_phase("3", "DONE")
}, error = function(e) {
  logmsg(sprintf("PHASE 3 ERROR: %s", conditionMessage(e)))
  mark_phase("3", paste0("FAILED: ", conditionMessage(e)))
})

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 4: Immune subcluster comprehensive case study
# ─────────────────────────────────────────────────────────────────────────────
logmsg("=== PHASE 4: Immune subcluster case study ===")
immune_case_csv <- file.path(OUT_ROOT, "immune_subcluster_casestudy.csv")

tryCatch({
  # Build cross-set gene → module assignment table (using large_wgcna as primary GGM)
  ggm_lw  <- read.csv(GGM_SETS[["GGM_large_wgcna"]], stringsAsFactors = FALSE)
  pb_wgcna <- read.csv(PB_SETS[["PB_wgcna"]], stringsAsFactors = FALSE)

  gene_to_ggm_lw  <- setNames(ggm_lw$top_module,  ggm_lw$gene_id)
  gene_to_pb_wgcna <- setNames(pb_wgcna$module,    pb_wgcna$gene_id)

  # For each immune subcluster, analyse marker gene concentration in modules
  immune_clusters <- sort(unique(mk_immune_raw$cluster))
  logmsg(sprintf("  Immune subclusters: %s", paste(immune_clusters, collapse = ",")))

  case_rows <- list()
  primer_genes <- NULL  # subcluster 4

  for (cl in immune_clusters) {
    cl_label <- paste0("ImmuneSubclust_c", cl)
    mk_genes_raw <- mk_immune_raw$gene[mk_immune_raw$cluster == cl]
    mk_atids <- unique(sapply(mk_genes_raw, resolve_to_atid, USE.NAMES = FALSE))
    mk_atids <- mk_atids[!is.na(mk_atids)]

    if (cl == 4) primer_genes <- mk_atids

    # 4a: Module assignment distribution
    ggm_assigns  <- gene_to_ggm_lw[mk_atids]
    ggm_assigns  <- ggm_assigns[!is.na(ggm_assigns)]
    pb_assigns   <- gene_to_pb_wgcna[mk_atids]
    pb_assigns   <- pb_assigns[!is.na(pb_assigns)]

    ggm_entropy <- if (length(ggm_assigns) > 0) {
      ggm_freq <- table(ggm_assigns) / length(ggm_assigns)
      -sum(ggm_freq * log(ggm_freq + 1e-12))
    } else NA_real_

    pb_entropy  <- if (length(pb_assigns) > 0) {
      pb_freq  <- table(pb_assigns)  / length(pb_assigns)
      -sum(pb_freq  * log(pb_freq  + 1e-12))
    } else NA_real_

    ggm_top_mod <- if (length(ggm_assigns) > 0) as.integer(names(which.max(table(ggm_assigns)))) else NA_integer_
    pb_top_mod  <- if (length(pb_assigns) > 0)  as.integer(names(which.max(table(pb_assigns))))  else NA_integer_
    ggm_top_frac <- if (length(ggm_assigns) > 0) max(table(ggm_assigns)) / length(ggm_assigns) else NA_real_
    pb_top_frac  <- if (length(pb_assigns) > 0)  max(table(pb_assigns))  / length(pb_assigns)  else NA_real_

    # 4b: Co-detection of marker group vs co-expression module
    mk_rna <- resolve_grp(mk_atids)
    sc_mk <- NULL
    if (length(mk_rna) >= 2L) {
      sc_mk <- tryCatch(
        score_gene_group(counts_mat, lognorm_mat, det_rate_all,
                         mk_rna, group_id = paste0(cl_label, "_markers"),
                         n_pair_cap = 2000L, seed = 98L, timeout_secs = 300L,
                         log_fn = logmsg),
        error = function(e) { logmsg(sprintf("  [c%d markers] %s", cl, conditionMessage(e))); NULL }
      )
    }

    case_rows[[as.character(cl)]] <- data.frame(
      immune_cluster  = cl,
      is_PRIMER       = (cl == 4),
      n_markers       = length(mk_atids),
      n_in_ggm_lw     = length(ggm_assigns),
      ggm_lw_top_module = ggm_top_mod,
      ggm_lw_top_frac = ggm_top_frac,
      ggm_lw_entropy  = ggm_entropy,
      n_in_pb_wgcna   = length(pb_assigns),
      pb_wgcna_top_module = pb_top_mod,
      pb_wgcna_top_frac = pb_top_frac,
      pb_wgcna_entropy = pb_entropy,
      marker_codet_gap = if (!is.null(sc_mk)) sc_mk$gap_codet else NA_real_,
      marker_codet_es  = if (!is.null(sc_mk)) sc_mk$effect_size_codet else NA_real_,
      stringsAsFactors = FALSE
    )

    logmsg(sprintf("  c%d: %d markers → GGM top_mod=%s (%.1f%%), PB top_mod=%s (%.1f%%)",
      cl, length(mk_atids),
      ifelse(is.na(ggm_top_mod), "NA", as.character(ggm_top_mod)),
      ifelse(is.na(ggm_top_frac), 0, ggm_top_frac * 100),
      ifelse(is.na(pb_top_mod), "NA", as.character(pb_top_mod)),
      ifelse(is.na(pb_top_frac), 0, pb_top_frac * 100)))
  }

  case_df <- do.call(rbind, case_rows)
  write.csv(case_df, immune_case_csv, row.names = FALSE, quote = TRUE)
  logmsg(sprintf("  Saved immune case study: %d rows", nrow(case_df)))

  # 4c: PRIMER cell (subcluster 4) full readout
  if (!is.null(primer_genes)) {
    primer_rows <- lapply(primer_genes, function(atid) {
      sym  <- id2sym[atid]
      ggm  <- gene_to_ggm_lw[atid]
      pb   <- gene_to_pb_wgcna[atid]
      ggm_kme <- if (!is.na(ggm)) {
        ggm_lw$kME[match(atid, ggm_lw$gene_id)]
      } else NA_real_
      pb_kme  <- if (!is.na(pb)) {
        pb_wgcna$kME[match(atid, pb_wgcna$gene_id)]
      } else NA_real_
      data.frame(
        gene_id = atid, display_label = make_label(atid),
        ggm_large_wgcna_module = ggm, ggm_large_wgcna_kME = ggm_kme,
        pb_wgcna_module = pb, pb_wgcna_kME = pb_kme,
        det_rate = det_rate_all[atid_to_rna[atid]],
        stringsAsFactors = FALSE
      )
    })
    primer_df <- do.call(rbind, primer_rows[!sapply(primer_rows, is.null)])
    primer_df <- primer_df[order(-primer_df$det_rate), ]
    write.csv(primer_df, file.path(OUT_ROOT, "PRIMER_readout.csv"),
              row.names = FALSE, quote = TRUE)
    logmsg(sprintf("  PRIMER readout: %d genes", nrow(primer_df)))
    logmsg(sprintf("  PRIMER GGM modules: %s", paste(sort(unique(na.omit(primer_df$ggm_large_wgcna_module))), collapse = ",")))
    logmsg(sprintf("  PRIMER PB modules:  %s", paste(sort(unique(na.omit(primer_df$pb_wgcna_module))), collapse = ",")))
  }

  # Module concentration figure for immune subclusters
  p_conc <- ggplot(case_df, aes(x = factor(immune_cluster), y = ggm_lw_top_frac,
                                  fill = factor(immune_cluster))) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_point(aes(y = pb_wgcna_top_frac), shape = 21, size = 2, fill = "white") +
    labs(x = "Immune subcluster", y = "Fraction of markers in top module",
         title = "Marker concentration in top co-expression module\n(bars=GGM; circles=PB)") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(FIG_DIR, "fig_immune_marker_concentration.png"), p_conc,
         width = 7, height = 4, dpi = 300, bg = "white")

  mark_phase("4", "DONE")
}, error = function(e) {
  logmsg(sprintf("PHASE 4 ERROR: %s", conditionMessage(e)))
  mark_phase("4", paste0("FAILED: ", conditionMessage(e)))
})

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 5: Single-cell visualization (UMAP feature plots)
# ─────────────────────────────────────────────────────────────────────────────
logmsg("=== PHASE 5: UMAP feature plots ===")
PURPLE_COLS <- c("lightgray", "#BFD3E6", "#9EBCDA", "#8C96C6",
                  "#8C6BB1", "#88419D", "#810F7C", "#4D004B")
CONDITIONS  <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")

tryCatch({
  # Condition column derivation (reuse featureplot_modules.R approach)
  sample2 <- seurat_obj$sample
  sample2 <- gsub("_rep[12]$", "", sample2)
  sample2 <- gsub("_(04|06|09|24)h$", "", sample2)
  sample2 <- gsub("^00_Mock$", "Mock", sample2)
  seurat_obj$condition <- factor(sample2, levels = CONDITIONS)

  umap_df <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$condition <- seurat_obj$condition

  plot_module_score_umap <- function(gene_atids, label, out_path) {
    rna_feats <- na.omit(atid_to_rna[gene_atids])
    if (length(rna_feats) < 2L) {
      logmsg(sprintf("  SKIP %s: < 2 resolvable genes", label))
      return(invisible(NULL))
    }
    tmp <- AddModuleScore(seurat_obj, features = list(rna_feats), name = "TmpMS_")
    scores <- tmp@meta.data[["TmpMS_1"]]
    rm(tmp)
    df <- umap_df
    df$score <- scores
    panels <- lapply(seq_along(CONDITIONS), function(k) {
      cond <- CONDITIONS[k]
      d    <- df[!is.na(df$condition) & df$condition == cond, ]
      d    <- d[order(d$score, na.last = FALSE), ]
      ttl  <- if (k == 1L) paste0(label, "\n(", cond, ")") else cond
      ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = score)) +
        geom_point(size = 0.4) +
        scale_color_gradientn(colours = PURPLE_COLS, na.value = "lightgray", name = "score") +
        labs(title = ttl) + NoAxes() +
        theme(plot.background = element_blank(), panel.background = element_blank(),
              plot.title = element_text(size = 8, face = "bold"))
    })
    fig <- wrap_plots(panels, nrow = 1)
    ggsave(out_path, fig, width = length(CONDITIONS) * 4.5, height = 4.5, dpi = 300, bg = "white")
    logmsg(sprintf("  Saved: %s", basename(out_path)))
    invisible(out_path)
  }

  # 5a: Co-expression modules most relevant to immune subclusters
  case_df2 <- read.csv(immune_case_csv, stringsAsFactors = FALSE)
  top_ggm_mods <- unique(na.omit(case_df2$ggm_lw_top_module))
  top_pb_mods  <- unique(na.omit(case_df2$pb_wgcna_top_module))

  for (gm in top_ggm_mods) {
    p_out <- file.path(FIG_DIR, sprintf("fig_GGM_lw_m%s_immune.png", gm))
    if (!file.exists(p_out)) {
      mod_genes <- ggm_lw$gene_id[ggm_lw$top_module == gm]
      tryCatch(plot_module_score_umap(mod_genes, paste0("GGM_lw_M", gm), p_out),
               error = function(e) logmsg(sprintf("  ERROR GGM m%s: %s", gm, conditionMessage(e))))
    }
  }

  # 5b: PRIMER module(s)
  if (!is.null(primer_genes)) {
    primer_ggm_mods <- unique(na.omit(gene_to_ggm_lw[primer_genes]))
    for (gm in primer_ggm_mods) {
      p_out <- file.path(FIG_DIR, sprintf("fig_PRIMER_GGM_m%s.png", gm))
      if (!file.exists(p_out)) {
        mod_genes <- ggm_lw$gene_id[ggm_lw$top_module == gm]
        tryCatch(plot_module_score_umap(mod_genes, paste0("PRIMER_GGM_M", gm), p_out),
                 error = function(e) logmsg(sprintf("  ERROR PRIMER GGM m%s: %s", gm, conditionMessage(e))))
      }
    }
    primer_rna <- na.omit(atid_to_rna[primer_genes])
    p_out <- file.path(FIG_DIR, "fig_PRIMER_marker_score.png")
    if (!file.exists(p_out))
      tryCatch(plot_module_score_umap(primer_genes, "PRIMER_markers", p_out),
               error = function(e) logmsg(sprintf("  ERROR PRIMER markers: %s", conditionMessage(e))))
  }

  # 5c: Contrasting examples (core vs mode-specific)
  core_rna <- resolve_grp(head(core_genes, 200))
  p_out <- file.path(FIG_DIR, "fig_CORE_score_example.png")
  if (!file.exists(p_out))
    tryCatch(plot_module_score_umap(head(core_genes, 200), "Core_genes_sample200", p_out),
             error = function(e) logmsg(sprintf("  ERROR core sample: %s", conditionMessage(e))))

  pb_only_rna <- resolve_grp(head(pb_only, 200))
  p_out <- file.path(FIG_DIR, "fig_PBONLY_score_example.png")
  if (!file.exists(p_out))
    tryCatch(plot_module_score_umap(head(pb_only, 200), "PB_only_sample200", p_out),
             error = function(e) logmsg(sprintf("  ERROR pb_only sample: %s", conditionMessage(e))))

  # 5d: Marker group vs matched co-expression module (2-3 immune subclusters)
  for (cl in head(immune_clusters, 3L)) {
    cl_label <- paste0("ImmuneSubclust_c", cl)
    mk_atids <- unique(sapply(mk_immune_raw$gene[mk_immune_raw$cluster == cl],
                               resolve_to_atid, USE.NAMES = FALSE))
    mk_atids <- mk_atids[!is.na(mk_atids)]

    # Marker score
    p_mk  <- file.path(FIG_DIR, sprintf("fig_immune_c%d_markers_score.png", cl))
    if (!file.exists(p_mk))
      tryCatch(plot_module_score_umap(mk_atids, paste0("c", cl, "_markers"), p_mk),
               error = function(e) logmsg(sprintf("  ERROR immune c%d markers: %s", cl, conditionMessage(e))))

    # Matched GGM module score
    ggm_top <- case_df2$ggm_lw_top_module[case_df2$immune_cluster == cl]
    if (length(ggm_top) > 0 && !is.na(ggm_top[1])) {
      mod_genes <- ggm_lw$gene_id[ggm_lw$top_module == ggm_top[1]]
      p_mod <- file.path(FIG_DIR, sprintf("fig_immune_c%d_GGM_m%s_score.png", cl, ggm_top[1]))
      if (!file.exists(p_mod))
        tryCatch(plot_module_score_umap(mod_genes, paste0("c", cl, "_GGM_M", ggm_top[1]), p_mod),
                 error = function(e) logmsg(sprintf("  ERROR immune c%d GGM mod: %s", cl, conditionMessage(e))))
    }
  }

  # FIGURES_MANIFEST.md
  figs <- list.files(FIG_DIR, pattern = "\\.png$", full.names = FALSE)
  man_lines <- c("# FIGURES MANIFEST", "", sprintf("Generated: %s", format(Sys.time())), "")
  for (f in sort(figs)) man_lines <- c(man_lines, paste0("- ", f))
  writeLines(man_lines, file.path(OUT_ROOT, "FIGURES_MANIFEST.md"))

  mark_phase("5", "DONE")
}, error = function(e) {
  logmsg(sprintf("PHASE 5 ERROR: %s", conditionMessage(e)))
  mark_phase("5", paste0("FAILED: ", conditionMessage(e)))
})

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 6: Integrated HTML report
# ─────────────────────────────────────────────────────────────────────────────
logmsg("=== PHASE 6: HTML report ===")
report_path <- file.path(OUT_ROOT, "CORE_MARKER_REPORT.html")

tryCatch({
  ensure_pkg("jsonlite")

  # Load saved CSVs (robust to which phases completed)
  bench_df2   <- if (file.exists(bench_csv)) read.csv(bench_csv, stringsAsFactors = FALSE) else data.frame()
  case_df3    <- if (file.exists(immune_case_csv)) read.csv(immune_case_csv, stringsAsFactors = FALSE) else data.frame()
  mult_df2    <- if (file.exists(file.path(OUT_ROOT, "marker_gene_multiplicity.csv")))
    read.csv(file.path(OUT_ROOT, "marker_gene_multiplicity.csv"), stringsAsFactors = FALSE) else data.frame()
  primer_df2  <- if (file.exists(file.path(OUT_ROOT, "PRIMER_readout.csv")))
    read.csv(file.path(OUT_ROOT, "PRIMER_readout.csv"), stringsAsFactors = FALSE) else data.frame()
  bench_sum2  <- if (file.exists(file.path(OUT_ROOT, "benchmark_summary.csv")))
    read.csv(file.path(OUT_ROOT, "benchmark_summary.csv"), stringsAsFactors = FALSE) else data.frame()
  overlap_df2 <- if (file.exists(file.path(OUT_ROOT, "core_marker_overlap.csv")))
    read.csv(file.path(OUT_ROOT, "core_marker_overlap.csv"), stringsAsFactors = FALSE) else data.frame()

  # Helper: embed PNG as base64
  embed_png <- function(path, width = "100%") {
    if (!file.exists(path)) return(sprintf("<p><em>[Figure not found: %s]</em></p>", basename(path)))
    b64 <- base64enc::dataURI(file = path, mime = "image/png")
    sprintf('<img src="%s" style="width:%s;max-width:900px;" />', b64, width)
  }

  # Helper: data.frame → HTML table
  df_to_html <- function(df, max_rows = 50L) {
    if (nrow(df) == 0) return("<p><em>No data.</em></p>")
    df <- head(df, max_rows)
    rows <- apply(df, 1, function(row) {
      paste0("<tr>", paste0("<td>", row, "</td>", collapse = ""), "</tr>")
    })
    paste0('<table class="tbl"><thead><tr>',
           paste0('<th>', colnames(df), '</th>', collapse = ''),
           '</tr></thead><tbody>',
           paste(rows, collapse = ''),
           '</tbody></table>')
  }

  # Plotly box trace JSON for benchmark
  make_plotly_box <- function(df, x_col, y_col, title, ylab) {
    if (nrow(df) == 0) return('<p><em>No data for plot.</em></p>')
    srcs <- unique(df[[x_col]])
    traces_json <- paste(sapply(srcs, function(s) {
      vals <- df[[y_col]][df[[x_col]] == s]
      sprintf('{"type":"box","name":%s,"y":[%s],"boxpoints":"outliers","marker":{"size":4}}',
        jsonlite::toJSON(s, auto_unbox=TRUE),
        paste(round(vals, 6), collapse = ","))
    }), collapse = ",")
    div_id <- paste0("plt_", gsub("[^a-zA-Z0-9]", "_", y_col))
    sprintf(
      '<div id="%s" style="width:800px;height:450px;"></div>
<script>Plotly.newPlot("%s",[%s],{title:{text:%s},yaxis:{title:%s},boxmode:"group"},
{displayModeBar:false,responsive:true});</script>',
      div_id, div_id, traces_json,
      jsonlite::toJSON(title, auto_unbox=TRUE),
      jsonlite::toJSON(ylab, auto_unbox=TRUE))
  }

  figs_available <- list.files(FIG_DIR, pattern = "\\.png$", full.names = TRUE)

  # PRIMER headline
  primer_headline <- if (nrow(primer_df2) > 0) {
    ggm_mods <- unique(na.omit(primer_df2$ggm_large_wgcna_module))
    pb_mods  <- unique(na.omit(primer_df2$pb_wgcna_module))
    sprintf("PRIMER cell (%d marker genes): best GGM module(s) = %s; PB module(s) = %s",
      nrow(primer_df2),
      paste(ggm_mods, collapse = ","),
      paste(pb_mods, collapse = ","))
  } else "PRIMER readout not available."

  # Benchmark headline
  bench_headline <- if (nrow(bench_sum2) > 0) {
    coexpr_row <- bench_sum2[grepl("Co-expression", bench_sum2$source_label, fixed=TRUE), ]
    marker_row  <- bench_sum2[grepl("Major", bench_sum2$source_label, fixed=TRUE), ]
    if (nrow(coexpr_row) > 0 && nrow(marker_row) > 0)
      sprintf("Co-expression modules: mean co-detection gap = %.4f (ES=%.2f) vs major markers: %.4f (ES=%.2f)",
        coexpr_row$mean_gap_codet[1], coexpr_row$mean_es_codet[1],
        marker_row$mean_gap_codet[1], marker_row$mean_es_codet[1])
    else "See benchmark_summary.csv"
  } else "Benchmark not completed."

  html_out <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Core/Marker Co-detection Report</title>
<script src="https://cdn.plot.ly/plotly-2.26.0.min.js" charset="utf-8"></script>
<style>
body{font-family:sans-serif;max-width:1100px;margin:auto;padding:1em;color:#222;}
h1{color:#1a3a5c;} h2{color:#2c5282;border-bottom:2px solid #bee3f8;padding-bottom:4px;}
h3{color:#2b6cb0;}
.tbl{border-collapse:collapse;width:100%;font-size:0.85em;}
.tbl th{background:#2c5282;color:#fff;padding:6px 8px;text-align:left;}
.tbl td{border:1px solid #ddd;padding:5px 8px;}
.tbl tr:nth-child(even){background:#f0f4ff;}
.note{background:#fffbe6;border-left:4px solid #f6ad55;padding:0.8em;margin:1em 0;}
.fig-wrap{text-align:center;margin:1em 0;}
section{margin-bottom:2.5em;}
</style>
</head>
<body>
<h1>Single-cell Co-detection Report: Marker Groups vs Co-expression Modules</h1>
<p><strong>Generated:</strong> ', format(Sys.time()), '</p>
<p><strong>Seurat object:</strong> combined_filtered.rds (', n_cells, ' cells, ', n_genes, ' genes)</p>

<section>
<h2>0. Executive summary</h2>
<p><strong>Benchmark:</strong> ', bench_headline, '</p>
<p><strong>PRIMER cell:</strong> ', primer_headline, '</p>
<p><strong>Core genes (1,095):</strong> Co-assigned by both GGM and pseudobulk modes —
see Phase 3 for expression character and marker overlap.</p>
<div class="note"><strong>Method:</strong> All groups scored on the same single-cell matrix.
Co-detection gap = within-group co-detection rate minus expression-frequency-matched-null
co-detection rate. Effect size = gap / SD(null). This controls for the confound that
highly-expressed genes co-detect at high rates by chance.</div>
</section>

<section>
<h2>1. The problem with marker gene lists</h2>
<p>Conventional snRNA analysis treats per-cluster marker lists as de-facto
co-expression groups. Three structural weaknesses:</p>
<ol>
<li><strong>Marker multiplicity</strong>: the same gene can be a significant marker for
multiple clusters (it is not an exclusive co-expression unit).</li>
<li><strong>Incomplete coverage</strong>: only genes with cluster-specific expression are
markers; broadly-expressed co-expressed genes are missed.</li>
<li><strong>Granularity dependence</strong>: marker groups change substantially with
clustering resolution.</li>
</ol>
<p>The co-detection benchmark below tests whether co-expression modules (which learn
gene groups directly from expression) overcome these weaknesses.</p>

<h3>Marker multiplicity (Phase 0)</h3>
', df_to_html(if(nrow(mult_df2)>0) mult_df2[order(-mult_df2$n_clusters), c("display_label","n_clusters")][seq_len(min(20,nrow(mult_df2))),] else data.frame()), '
<div class="fig-wrap">', embed_png(file.path(FIG_DIR, "fig_marker_multiplicity.png"), "500px"), '</div>
</section>

<section>
<h2>2. The co-detection benchmark (Phase 2)</h2>
<div class="note"><strong>Design:</strong> Three sources evaluated on identical metrics:
(i) major-cluster marker groups, (ii) immune-subcluster marker groups,
(iii) co-expression modules (all 6 sets). Same single-cell matrix. Matched-null controls
for detection-rate confound.</div>

<h3>Summary by source</h3>
', df_to_html(bench_sum2), '
', if(nrow(bench_df2)>0) make_plotly_box(bench_df2[!is.na(bench_df2$gap_codet),],
     "source", "gap_codet",
     "Co-detection gap by source", "Co-detection gap (within − null)") else '', '

<h3>Granularity: major vs immune-subcluster markers</h3>
<div class="fig-wrap">', embed_png(file.path(FIG_DIR, "fig_granularity_effect.png"), "500px"), '</div>

<div class="note"><strong>Confounds:</strong> Co-expression modules are typically larger than
marker groups. Effect sizes (gap / SD_null) are the appropriate comparison, not raw gaps.
The matched-null already controls for detection-rate differences between genes.</div>
</section>

<section>
<h2>3. Core vs mode-specific genes (Phase 3)</h2>
<p>1,095 "core" genes co-assigned by both GGM and pseudobulk modes; 8,371 GGM-only;
4,048 PB-only (from downstream/core_vs_modespecific.csv).</p>

<h3>Marker overlap (hypergeometric)</h3>
', df_to_html(overlap_df2), '

<h3>Expression character</h3>
<div class="fig-wrap">', embed_png(file.path(FIG_DIR, "fig_core_vs_modespecific_expression.png"), "500px"), '</div>
<div class="fig-wrap">', embed_png(file.path(FIG_DIR, "fig_core_marker_overlap.png"), "500px"), '</div>
<div class="fig-wrap">', embed_png(file.path(FIG_DIR, "fig_core_codetection.png"), "500px"), '</div>
</section>

<section>
<h2>4. Immune subcluster case study (Phase 4)</h2>
<p>Immune-active subclusters (major 3/7/11 pooled, subclustered). PRIMER cell = subcluster 4.</p>

<h3>Module concentration per subcluster</h3>
<div class="note">High "top_frac" = marker genes concentrate in one co-expression module
(coherent). Low = markers scatter across modules (incoherent with co-expression structure).</div>
', df_to_html(case_df3), '
<div class="fig-wrap">', embed_png(file.path(FIG_DIR, "fig_immune_marker_concentration.png"), "700px"), '</div>

<h3>PRIMER cell (subcluster 4)</h3>
<p>', primer_headline, '</p>
', df_to_html(head(primer_df2, 30)), '
</section>

<section>
<h2>5. Single-cell visualization (Phase 5)</h2>
',
paste(sapply(sort(figs_available), function(fp) {
  nm <- basename(fp)
  paste0('<div class="fig-wrap"><p><strong>', nm, '</strong></p>', embed_png(fp, "900px"), '</div>')
}), collapse = "\n"),
'
</section>

<section>
<h2>6. Summary and limitations</h2>
<p>', bench_headline, '</p>
<ul>
<li>Co-expression modules learn groups directly from expression; markers from
cluster-differential expression. The two are complementary, not identical.</li>
<li>Matched-null controls for detection-rate confound but NOT for gene-count
differences between groups (use effect sizes for cross-source comparison).</li>
<li>PRIMER cell: its marker genes are analysed for module assignment concentration;
the co-expression module capturing the PRIMER signature is the one with highest
marker concentration.</li>
<li>Phase 3 data show whether core genes are more broadly expressed
(high detection rate) vs mode-specific genes (context-restricted) — see figures.</li>
</ul>
<p><strong>Reusability:</strong> The co-detection engine
(inst/scripts/codetection_eval.R) is generic — accepts any gene group + sparse
matrix. No prior-knowledge dependency.</p>
</section>

</body>
</html>')

  writeLines(html_out, report_path)
  report_size <- file.size(report_path)
  logmsg(sprintf("  HTML report: %s (%.1f MB)", report_path, report_size / 1e6))

  mark_phase("6", "DONE")
}, error = function(e) {
  logmsg(sprintf("PHASE 6 ERROR: %s", conditionMessage(e)))
  mark_phase("6", paste0("FAILED: ", conditionMessage(e)))
})

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 7: Final console report
# ─────────────────────────────────────────────────────────────────────────────
logmsg("=== PHASE 7: Final report ===")
logmsg(sprintf("Seurat: %d cells, %d genes", n_cells, n_genes))
logmsg(sprintf("Marker files: %d major clusters, %d immune subclusters",
  length(marker_groups_major), length(marker_groups_immune)))
logmsg("Phase status:")
for (nm in names(phase_status))
  logmsg(sprintf("  PHASE %s: %s", nm, phase_status[[nm]]))

n_figs <- length(list.files(FIG_DIR, pattern = "\\.png$"))
logmsg(sprintf("Figures saved: %d → %s", n_figs, FIG_DIR))
logmsg(sprintf("Report: %s", report_path))
logmsg(sprintf("Benchmark CSV: %s", bench_csv))
logmsg(sprintf("PRIMER readout: %s", file.path(OUT_ROOT, "PRIMER_readout.csv")))
logmsg("=== core_marker_analysis.R complete ===")

close(log_con)
