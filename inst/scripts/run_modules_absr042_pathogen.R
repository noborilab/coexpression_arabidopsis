## Pseudobulk co-expression modules at mean |r| >= 0.42 (pathogen_multiome)
## Uses existing robustness_result.rds from pseudobulk_zscore_spearman pipeline.
## Louvain = PRIMARY; WGCNA power=1 = SECONDARY (20-min hard timeout).
##
## Steps:
##   1. Network stats at |r| >= 0.42 (confirm ~752k pairs / ~5,450 genes)
##   2. Louvain modules → modules_absr042/large_louvain/
##   3. WGCNA quick test → modules_absr042/large_wgcna_test/ (or timeout record)
##   4. Feature plots (Louvain) + BON3/WRKY readout

suppressPackageStartupMessages({
  # Load from source if running from repo root (picks up R/interpret.R changes)
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(CoexprArabidopsis)
  }
  library(igraph)
  library(ggplot2)
  library(patchwork)
})

t_global <- proc.time()

# ==============================================================================
# PARAMETERS
# ==============================================================================

MIN_ABS_R  <- 0.42
CONDITIONS <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")
DATASET_ID <- "pathogen_multiome"

SEURAT_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/",
  "Projects/SA_PTI_ETI_single_cell/",
  "SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed/",
  "combined_filtered.rds"
)
TF_META_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Nobori Lab (TSL) Team Folder/",
  "shared/datasets/from_Ben/for_tatsuya/data/motifs-2026/",
  "Athaliana_motifs_metadata.tsv"
)

RESULTS_DIR     <- file.path("results", DATASET_ID)
OUT_DIR         <- file.path(RESULTS_DIR, "pseudobulk_zscore_spearman")
OUT_MOD_DIR     <- file.path(OUT_DIR, "modules_absr042")
PLOTS_DIR       <- file.path(OUT_MOD_DIR, "plots")
SYMBOL_MAP_PATH <- file.path(RESULTS_DIR, "symbol_map.csv")
WRKY_PATH       <- file.path(RESULTS_DIR, "geneset_lookups", "WRKY_GGM_vs_PB.csv")
BON3_ID         <- "AT1G08860"

dir.create(OUT_MOD_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Load shared resources
# ==============================================================================

message("\n=== Loading shared resources ===")

symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)
sym_lookup  <- setNames(symbol_map$gene_symbol, symbol_map$gene_id)
message("symbol_map: ", nrow(symbol_map), " entries")

rob <- readRDS(file.path(OUT_DIR, "robustness_result.rds"))
ps  <- rob$pair_scores
message("pair_scores loaded: ", nrow(ps), " pairs")

.add_symbol <- function(df, col = "gene_id") {
  if (col %in% names(df)) df$gene_symbol <- sym_lookup[df[[col]]]
  df
}

# Per-condition edge proxy from I_ columns (for annotate_context)
# Uses z_bar as per-condition weight proxy — sufficient to rank conditions.
message("Building per-condition edge proxy from I_ columns...")
network_list_ctx <- setNames(lapply(CONDITIONS, function(cond) {
  i_col <- paste0("I_", cond)
  if (!i_col %in% names(ps))
    return(list(edge_table = data.frame(gene_id_A = character(), gene_id_B = character(),
                                        weight = numeric(), stringsAsFactors = FALSE)))
  ps_c <- ps[!is.na(ps[[i_col]]) & ps[[i_col]] == 1L, , drop = FALSE]
  list(edge_table = data.frame(
    gene_id_A = ps_c$gene_id_A,
    gene_id_B = ps_c$gene_id_B,
    weight    = abs(tanh(ps_c$z_bar)),
    stringsAsFactors = FALSE
  ))
}), CONDITIONS)

message("network_list_ctx sizes: ",
        paste(vapply(CONDITIONS, function(c)
          sprintf("%s=%d", c, nrow(network_list_ctx[[c]]$edge_table)),
          character(1)), collapse = " | "))

# ==============================================================================
# Louvain module builder (min_abs_r filter)
# ==============================================================================

.build_louvain_absr <- function(rob, min_abs_r, min_module_size = 30L) {
  ps <- rob$pair_scores
  ps <- ps[!is.na(ps$z_bar) & abs(tanh(ps$z_bar)) >= min_abs_r, , drop = FALSE]
  if (nrow(ps) == 0L)
    stop("No pairs after min_abs_r filter = ", min_abs_r)

  edges <- data.frame(
    from   = ps$gene_id_A,
    to     = ps$gene_id_B,
    weight = pmin(abs(tanh(ps$z_bar)), 1.0),
    stringsAsFactors = FALSE
  )

  g  <- igraph::graph_from_data_frame(edges, directed = FALSE)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  memb     <- igraph::membership(cl)
  gene_ids <- names(memb)
  top_lbl  <- as.integer(memb)

  comm_sz <- table(top_lbl)
  small   <- as.integer(names(comm_sz[comm_sz < min_module_size]))
  top_lbl[top_lbl %in% small] <- 0L

  live    <- sort(unique(top_lbl[top_lbl > 0L]))
  relabel <- setNames(seq_along(live), as.character(live))
  top_lbl <- ifelse(top_lbl == 0L, 0L,
                    as.integer(relabel[as.character(top_lbl)]))

  all_genes <- sort(union(edges$from, edges$to))
  n         <- length(all_genes)
  gene_idx  <- setNames(seq_along(all_genes), all_genes)

  A <- matrix(0.0, n, n, dimnames = list(all_genes, all_genes))
  iA <- gene_idx[edges$from]; iB <- gene_idx[edges$to]
  A[cbind(iA, iB)] <- edges$weight
  A[cbind(iB, iA)] <- edges$weight

  tl_ord      <- setNames(top_lbl, gene_ids)[all_genes]
  unique_mods <- sort(unique(tl_ord[tl_ord > 0L]))

  kme_vec <- rep(NA_real_, n)
  for (m in unique_mods) {
    mod_idx  <- which(tl_ord == m)
    if (length(mod_idx) < 2L) next
    mod_mean <- colMeans(A[mod_idx, , drop = FALSE])
    for (i in mod_idx) kme_vec[i] <- cor(A[i, ], mod_mean)
  }

  gene_module <- data.frame(
    gene_id    = all_genes,
    top_module = as.integer(tl_ord),
    sub_module = NA_integer_,
    kME        = kme_vec,
    stringsAsFactors = FALSE
  )

  mod_counts <- table(gene_module$top_module[gene_module$top_module > 0L])
  module_meta <- data.frame(
    module_id              = as.integer(names(mod_counts)),
    n_genes                = as.integer(mod_counts),
    label                  = NA_character_,
    top_organ_or_condition = NA_character_,
    delta_treatment        = NA_character_,
    go_top                 = NA_character_,
    zsummary               = NA_real_,
    preservation_method    = NA_character_,
    stringsAsFactors = FALSE
  )

  module_hier <- data.frame(sub_module = integer(), top_module = integer(),
                            stringsAsFactors = FALSE)

  hub_list <- lapply(as.integer(names(mod_counts)), function(m) {
    sg <- gene_module[gene_module$top_module == m & !is.na(gene_module$kME), ]
    sg <- head(sg[order(sg$kME, decreasing = TRUE), ], 20L)
    if (nrow(sg) == 0L) return(NULL)
    data.frame(module_id = m, gene_id = sg$gene_id, gene_symbol = NA_character_,
               kME = sg$kME, hub_rank = seq_len(nrow(sg)), stringsAsFactors = FALSE)
  })
  hub_genes <- do.call(rbind, Filter(Negate(is.null), hub_list))
  if (is.null(hub_genes) || nrow(hub_genes) == 0L)
    hub_genes <- data.frame(module_id = integer(), gene_id = character(),
                            gene_symbol = character(), kME = numeric(),
                            hub_rank = integer(), stringsAsFactors = FALSE)

  module_tfs <- data.frame(module_id = integer(), gene_id = character(),
                           gene_symbol = character(), tf_family = character(),
                           stringsAsFactors = FALSE)

  ME_mat <- matrix(NA_real_, nrow = length(unique_mods), ncol = n,
                   dimnames = list(paste0("ME", unique_mods), all_genes))
  for (m in unique_mods) {
    mi <- which(tl_ord == m)
    ME_mat[paste0("ME", m), ] <- colMeans(A[mi, , drop = FALSE])
  }
  MEs <- as.data.frame(t(ME_mat))

  list(gene_module = gene_module, module_meta = module_meta,
       module_hier = module_hier, hub_genes = hub_genes,
       module_tfs = module_tfs, eigengenes = MEs)
}

# Annotate and save one module set (preservation skipped)
.annotate_and_save_absr <- function(set_name, mod_input, network_list_ctx,
                                    tf_meta_path, out_mod_dir) {
  t0      <- proc.time()
  set_dir <- file.path(out_mod_dir, set_name)
  dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)
  message("  Annotating: ", set_name)

  mod_input$gene_module <- .add_symbol(mod_input$gene_module)
  mod_input$hub_genes   <- .add_symbol(mod_input$hub_genes)

  mod_input <- tryCatch(
    annotate_context(mod_input, network_list_ctx),
    error = function(e) { message("  annotate_context failed: ", conditionMessage(e)); mod_input }
  )

  mod_input <- tryCatch(
    annotate_go(mod_input, org_db = "org.At.tair.db", pval_cut = 0.05),
    error = function(e) { message("  annotate_go failed: ", conditionMessage(e)); mod_input }
  )

  if (file.exists(tf_meta_path)) {
    mod_input <- tryCatch(
      annotate_tfs(mod_input, tf_meta_path),
      error = function(e) { message("  annotate_tfs failed: ", conditionMessage(e)); mod_input }
    )
  } else {
    message("  TF file not found: ", tf_meta_path)
  }

  # preservation skipped (skip_preservation = TRUE)

  write.csv(mod_input$gene_module,  file.path(set_dir, "gene_module.csv"),  row.names = FALSE)
  write.csv(mod_input$module_meta,  file.path(set_dir, "module_meta.csv"),  row.names = FALSE)
  write.csv(mod_input$module_hier,  file.path(set_dir, "module_hier.csv"),  row.names = FALSE)
  write.csv(mod_input$hub_genes,    file.path(set_dir, "hub_genes.csv"),    row.names = FALSE)
  write.csv(mod_input$module_tfs,   file.path(set_dir, "module_tfs.csv"),   row.names = FALSE)
  write.csv(as.data.frame(mod_input$eigengenes),
            file.path(set_dir, "eigengenes.csv"), row.names = TRUE)
  saveRDS(mod_input, file.path(set_dir, "module_input.rds"))

  gm        <- mod_input$gene_module
  n_mods    <- length(unique(gm$top_module[gm$top_module > 0L]))
  n_grey    <- sum(is.na(gm$top_module) | gm$top_module == 0L)
  sizes     <- sort(as.integer(table(gm$top_module[gm$top_module > 0L])))
  size_min  <- if (length(sizes) > 0L) min(sizes)    else NA_integer_
  size_med  <- if (length(sizes) > 0L) median(sizes) else NA_real_
  size_max  <- if (length(sizes) > 0L) max(sizes)    else NA_integer_

  message(sprintf(
    "  [%s] modules=%d | grey=%d/%d (%.1f%%) | size min/med/max=%d/%.0f/%d | %.1f min",
    set_name, n_mods, n_grey, nrow(gm), 100 * n_grey / nrow(gm),
    size_min, size_med, size_max,
    (proc.time() - t0)[["elapsed"]] / 60
  ))

  mod_input
}

# ==============================================================================
# STEP 1: Network stats at |r| >= MIN_ABS_R
# ==============================================================================

message("\n==== Step 1: Network stats at mean |r| >= ", MIN_ABS_R, " ====")
t1 <- proc.time()

ps042   <- ps[!is.na(ps$z_bar) & abs(tanh(ps$z_bar)) >= MIN_ABS_R, , drop = FALSE]
n_pairs <- nrow(ps042)
n_genes <- length(unique(c(ps042$gene_id_A, ps042$gene_id_B)))
density <- n_pairs / (n_genes * (n_genes - 1L) / 2)

message(sprintf("  n_pairs = %d  |  n_genes = %d  |  density = %.4f",
                n_pairs, n_genes, density))

if (abs(n_pairs - 752000) / 752000 > 0.25)
  warning("n_pairs (", n_pairs, ") deviates >25% from expected ~752k")
if (abs(n_genes - 5450) / 5450 > 0.25)
  warning("n_genes (", n_genes, ") deviates >25% from expected ~5,450")

message(sprintf("Step 1 complete — n_pairs=%d  n_genes=%d  density=%.4f",
                n_pairs, n_genes, density))

# ==============================================================================
# STEP 2: Louvain modules (PRIMARY)
# ==============================================================================

message("\n==== Step 2: Louvain modules (primary) ====")
t2 <- proc.time()

n_mods_l <- NA_integer_; n_grey_l <- NA_integer_
sizes_l  <- integer(0);  louvain_result <- NULL

message("  Building Louvain graph at |r| >= ", MIN_ABS_R, "...")

louvain_result <- tryCatch(
  .build_louvain_absr(rob, min_abs_r = MIN_ABS_R, min_module_size = 30L),
  error = function(e) { message("  ERROR in Louvain build: ", conditionMessage(e)); NULL }
)

if (!is.null(louvain_result)) {
  gm_l      <- louvain_result$gene_module
  n_mods_l  <- length(unique(gm_l$top_module[gm_l$top_module > 0L]))
  n_grey_l  <- sum(is.na(gm_l$top_module) | gm_l$top_module == 0L)
  sizes_l   <- sort(as.integer(table(gm_l$top_module[gm_l$top_module > 0L])))
  message(sprintf("  Louvain raw: %d modules | grey=%d/%d (%.1f%%)",
                  n_mods_l, n_grey_l, nrow(gm_l), 100 * n_grey_l / nrow(gm_l)))
  message(sprintf("  Size: min=%d | median=%.0f | max=%d",
                  min(sizes_l), median(sizes_l), max(sizes_l)))

  louvain_result <- .annotate_and_save_absr(
    set_name         = "large_louvain",
    mod_input        = louvain_result,
    network_list_ctx = network_list_ctx,
    tf_meta_path     = TF_META_PATH,
    out_mod_dir      = OUT_MOD_DIR
  )
}

message(sprintf(
  "Step 2 complete — Louvain: modules=%s  grey_rate=%s  size_min/med/max=%s/%s/%s (%.1f min)",
  if (is.na(n_mods_l)) "FAILED" else as.character(n_mods_l),
  if (is.na(n_grey_l) || is.null(louvain_result)) "NA" else
    sprintf("%.1f%%", 100 * n_grey_l / nrow(louvain_result$gene_module)),
  if (length(sizes_l) == 0L) "NA" else as.character(min(sizes_l)),
  if (length(sizes_l) == 0L) "NA" else sprintf("%.0f", median(sizes_l)),
  if (length(sizes_l) == 0L) "NA" else as.character(max(sizes_l)),
  (proc.time() - t2)[["elapsed"]] / 60
))

# ==============================================================================
# STEP 3: WGCNA quick test (SECONDARY — 20-min hard timeout)
# ==============================================================================

message("\n==== Step 3: WGCNA quick test (secondary, 20-min timeout) ====")
t3 <- proc.time()

wgcna_result   <- NULL
wgcna_outcome  <- "not_run"
n_mods_w  <- NA_integer_; n_grey_w <- NA_integer_

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  wgcna_outcome <- "wgcna_not_installed"
  message("  WGCNA not installed; skipping WGCNA test.")
} else {
  message("  Launching WGCNA build (power=1, timeout=20 min)...")
  setTimeLimit(elapsed = 1200)
  wgcna_result <- tryCatch({
    r <- build_wgcna_modules(
      rob             = rob,
      network_list    = list(),   # not used inside build_wgcna_modules
      min_abs_r       = MIN_ABS_R,
      soft_power      = 1L,
      merge_cut       = 0.25,
      min_module_size = 30L,
      sub_merge_cut   = 0.10
    )
    setTimeLimit()   # clear limit on success
    r
  }, error = function(e) {
    setTimeLimit()   # always clear
    if (grepl("elapsed time limit", conditionMessage(e), ignore.case = TRUE)) {
      wgcna_outcome <<- "timeout"
      message("  WGCNA build exceeded 20 min timeout at |r|>=", MIN_ABS_R)
      writeLines(
        sprintf("WGCNA build exceeded 20 min timeout at |r|>=%s\n%s",
                MIN_ABS_R, format(Sys.time())),
        file.path(OUT_MOD_DIR, "large_wgcna_test_timeout.txt")
      )
    } else {
      wgcna_outcome <<- "error"
      message("  WGCNA build failed: ", conditionMessage(e))
    }
    NULL
  })
  setTimeLimit()  # belt-and-suspenders clear
}

if (!is.null(wgcna_result)) {
  gm_w     <- wgcna_result$gene_module
  n_mods_w <- length(unique(gm_w$top_module[gm_w$top_module > 0L]))
  n_grey_w <- sum(is.na(gm_w$top_module) | gm_w$top_module == 0L)
  sizes_w  <- sort(as.integer(table(gm_w$top_module[gm_w$top_module > 0L])))

  if (n_mods_w <= 1L) {
    wgcna_outcome <- "collapsed"
    message(sprintf(
      "  WGCNA DIAGNOSTIC: collapsed to %d module(s) — useful signal, not a failure",
      n_mods_w))
  } else {
    wgcna_outcome <- "completed"
  }

  message(sprintf("  WGCNA: %d modules | grey=%d/%d (%.1f%%) | size min/med/max=%d/%.0f/%d",
                  n_mods_w, n_grey_w, nrow(gm_w), 100 * n_grey_w / nrow(gm_w),
                  min(sizes_w), median(sizes_w), max(sizes_w)))

  wgcna_result <- .annotate_and_save_absr(
    set_name         = "large_wgcna_test",
    mod_input        = wgcna_result,
    network_list_ctx = network_list_ctx,
    tf_meta_path     = TF_META_PATH,
    out_mod_dir      = OUT_MOD_DIR
  )
}

message(sprintf("Step 3 complete — WGCNA: %s (%.1f min)", wgcna_outcome,
                (proc.time() - t3)[["elapsed"]] / 60))

# ==============================================================================
# STEP 4: Feature plots (Louvain) + BON3/WRKY readout
# ==============================================================================

message("\n==== Step 4: Feature plots + BON3/WRKY readout ====")
t4 <- proc.time()

PURPLE_COLS <- c("lightgray", "#BFD3E6", "#9EBCDA", "#8C96C6",
                 "#8C6BB1", "#88419D", "#810F7C", "#4D004B")
TOP_N <- 4L

# --- Feature plots ---
if (!is.null(louvain_result) && file.exists(SEURAT_PATH)) {
  message("  Loading Seurat object for feature plots...")
  suppressPackageStartupMessages(library(Seurat))

  sobj <- tryCatch(readRDS(SEURAT_PATH), error = function(e) {
    message("  Failed to load Seurat: ", conditionMessage(e)); NULL
  })

  if (!is.null(sobj)) {
    s2   <- sobj$sample
    s2   <- gsub("_rep[12]$", "", s2)
    s2   <- gsub("_(04|06|09|24)h$", "", s2)
    s2   <- gsub("^00_Mock$", "Mock", s2)
    sobj$condition <- factor(s2, levels = CONDITIONS)
    DefaultAssay(sobj) <- "RNA"
    RNA_GENES <- rownames(sobj[["RNA"]])

    .resolve_feature <- function(gene_id, gene_symbol = NA_character_) {
      cands <- c(
        if (!is.na(gene_symbol) && nchar(trimws(gene_symbol)) > 0)
          gene_symbol else character(0),
        gene_id,
        { sy <- symbol_map$gene_symbol[symbol_map$gene_id == gene_id]
          if (length(sy) > 0L && !is.na(sy[1L])) sy[1L] else character(0) }
      )
      for (cand in cands) if (cand %in% RNA_GENES) return(cand)
      NA_character_
    }

    .make_label <- function(gid, gsym)
      if (!is.na(gsym) && nchar(trimws(gsym)) > 0L) paste0(gid, " / ", gsym) else gid

    hub_df   <- louvain_result$hub_genes
    gm_df    <- louvain_result$gene_module
    all_mods <- sort(unique(hub_df$module_id))
    pdf_path <- file.path(PLOTS_DIR, "featureplots_large_louvain.pdf")

    message(sprintf("  large_louvain: %d modules → %s", length(all_mods), pdf_path))

    pdf(pdf_path, width = length(CONDITIONS) * 5, height = (TOP_N + 1L) * 5, onefile = TRUE)

    for (mid in all_mods) {
      top_genes <- hub_df[hub_df$module_id == mid, ]
      top_genes <- head(top_genes[order(top_genes$hub_rank), ], TOP_N)
      if (nrow(top_genes) == 0L) next

      top_genes$feature <- mapply(
        .resolve_feature, top_genes$gene_id,
        ifelse(is.na(top_genes$gene_symbol), NA_character_, top_genes$gene_symbol),
        SIMPLIFY = TRUE
      )
      top_genes_ok <- top_genes[!is.na(top_genes$feature), , drop = FALSE]
      if (nrow(top_genes_ok) == 0L) next

      all_ids <- gm_df$gene_id[gm_df$top_module == mid]
      all_rna <- vapply(all_ids, function(gid) {
        sy <- symbol_map$gene_symbol[symbol_map$gene_id == gid]
        .resolve_feature(gid, if (length(sy) > 0L && !is.na(sy[1L])) sy[1L] else NA_character_)
      }, character(1))
      all_rna <- na.omit(all_rna)

      tryCatch({
        tmp <- AddModuleScore(sobj, features = list(all_rna), name = "TmpMS_")
        module_scores <- tmp@meta.data[["TmpMS_1"]]
        rm(tmp)

        umap_df <- as.data.frame(Embeddings(sobj, reduction = "umap"))
        colnames(umap_df) <- c("UMAP_1", "UMAP_2")
        umap_df$condition    <- sobj$condition
        umap_df$module_score <- module_scores

        gene_rows <- lapply(seq_len(nrow(top_genes_ok)), function(i) {
          feat  <- top_genes_ok$feature[i]
          label <- .make_label(top_genes_ok$gene_id[i], top_genes_ok$gene_symbol[i])
          fp <- FeaturePlot(sobj, features = feat, split.by = "condition",
                            reduction = "umap", pt.size = 0.5, order = TRUE,
                            max.cutoff = "q99", min.cutoff = "q1")
          fp <- fp & scale_colour_gradientn(colours = PURPLE_COLS, na.value = "lightgray") &
            NoAxes()
          for (k in seq_along(fp)) {
            new_title <- if (k == 1L) paste0(label, "\n(", CONDITIONS[k], ")") else CONDITIONS[k]
            fp[[k]] <- fp[[k]] + labs(title = new_title) +
              theme(plot.title = element_text(size = 9, face = "bold"))
          }
          fp
        })

        score_panels <- lapply(seq_along(CONDITIONS), function(k) {
          cond_k <- CONDITIONS[k]
          df_k   <- umap_df[umap_df$condition == cond_k, ]
          df_k   <- df_k[order(df_k$module_score, na.last = FALSE), ]
          title  <- if (k == 1L) paste0("M", mid, " score\n(", cond_k, ")") else cond_k
          ggplot(df_k, aes(x = UMAP_1, y = UMAP_2, color = module_score)) +
            geom_point(size = 0.5) +
            scale_color_gradientn(colours = PURPLE_COLS, na.value = "lightgray", name = "score") +
            labs(title = title) + NoAxes() +
            theme(plot.background = element_blank(), panel.background = element_blank(),
                  plot.title = element_text(size = 9, face = "bold"))
        })

        print(wrap_plots(c(gene_rows, list(wrap_plots(score_panels, nrow = 1L))), ncol = 1L))

      }, error = function(e)
        message(sprintf("  ERROR M%s: %s", mid, conditionMessage(e)))
      )
    }

    dev.off()
    message(sprintf("  Saved: %s", pdf_path))
    rm(sobj); gc()
  }
} else if (is.null(louvain_result)) {
  message("  Louvain failed — skipping feature plots.")
} else {
  message("  Seurat object not found at expected path — skipping feature plots.")
  message("  Expected: ", SEURAT_PATH)
}

# --- BON3 / WRKY readout ---

message("\n  --- BON3 / WRKY readout ---")

bon3_all <- ps[ps$gene_id_A == BON3_ID | ps$gene_id_B == BON3_ID, , drop = FALSE]
bon3_all$partner <- ifelse(bon3_all$gene_id_A == BON3_ID, bon3_all$gene_id_B, bon3_all$gene_id_A)
bon3_042 <- bon3_all[!is.na(bon3_all$z_bar) & abs(tanh(bon3_all$z_bar)) >= MIN_ABS_R, ]
bon3_042  <- bon3_042[order(abs(tanh(bon3_042$z_bar)), decreasing = TRUE), ]

cat("\n=== BON3 / WRKY Summary (|r|>=", MIN_ABS_R, ", large_louvain) ===\n\n", sep = "")
cat(sprintf("BON3 (%s):\n", BON3_ID))
cat(sprintf("  Partners in pair_scores (any |r|): %d\n", nrow(bon3_all)))
cat(sprintf("  Partners at |r|>=%.2f: %d\n", MIN_ABS_R, nrow(bon3_042)))

if (!is.null(louvain_result)) {
  gm_l <- louvain_result$gene_module
  r_bon3 <- gm_l[gm_l$gene_id == BON3_ID, , drop = FALSE]
  if (nrow(r_bon3) == 0L) {
    cat(sprintf("  BON3: not in Louvain network (no edges at |r|>=%.2f)\n", MIN_ABS_R))
  } else {
    top_mod <- r_bon3$top_module[1L]
    mod_lbl <- if (!is.na(top_mod) && top_mod == 0L) "grey(0)" else as.character(top_mod)
    cat(sprintf("  BON3: module=%s | kME=%.5f | n_partners_at_042=%d\n",
                mod_lbl, r_bon3$kME[1L], nrow(bon3_042)))
  }
}

if (nrow(bon3_042) > 0L) {
  top10 <- head(bon3_042, 10L)
  top10$partner_symbol <- sym_lookup[top10$partner]
  cat(sprintf("\n  BON3 top-10 partners at |r|>=%.2f (by |r|):\n", MIN_ABS_R))
  for (i in seq_len(nrow(top10))) {
    sym_s <- ifelse(is.na(top10$partner_symbol[i]), "NA", top10$partner_symbol[i])
    cat(sprintf("    %2d. %s (%s)  |r|=%.4f  z_bar=%.4f  R_score=%.0f\n",
                i, top10$partner[i], sym_s,
                abs(tanh(top10$z_bar[i])), top10$z_bar[i], top10$R_score[i]))
  }
}

# WRKY family
wrky_ids <- character(0)
if (file.exists(WRKY_PATH)) {
  wrky_df  <- read.csv(WRKY_PATH, stringsAsFactors = FALSE)
  wrky_ids <- wrky_df$gene_id
  message("\n  WRKY IDs from WRKY_GGM_vs_PB.csv: ", length(wrky_ids))
} else if (file.exists(TF_META_PATH)) {
  tf_meta <- read.table(TF_META_PATH, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, quote = "")
  fam_col <- if      ("tf_family" %in% names(tf_meta)) "tf_family"
             else if ("family"    %in% names(tf_meta)) "family"
             else if ("class"     %in% names(tf_meta)) "class"
             else NA_character_
  if (!is.na(fam_col)) {
    wrky_ids <- tf_meta$motif_id[grepl("WRKY", tf_meta[[fam_col]], ignore.case = TRUE)]
    message("\n  WRKY IDs from TF metadata: ", length(wrky_ids))
  }
} else {
  message("\n  WARNING: no WRKY source found (WRKY_PATH and TF_META_PATH both missing)")
}

cat(sprintf("\nWRKY family (%d AT-IDs) in large_louvain:\n", length(wrky_ids)))

if (length(wrky_ids) > 0L && !is.null(louvain_result)) {
  gm_l     <- louvain_result$gene_module
  gm_wrky  <- gm_l[gm_l$gene_id %in% wrky_ids, , drop = FALSE]
  n_w_tot  <- nrow(gm_wrky)
  n_w_asn  <- sum(!is.na(gm_wrky$top_module) & gm_wrky$top_module > 0L)
  n_w_grey <- n_w_tot - n_w_asn

  kme_a <- gm_wrky$kME[!is.na(gm_wrky$kME) & !is.na(gm_wrky$top_module) &
                        gm_wrky$top_module > 0L]
  kme_str <- if (length(kme_a) > 0L)
    sprintf("kME: min=%.3f med=%.3f max=%.3f", min(kme_a), median(kme_a), max(kme_a))
  else "kME: N/A"

  cat(sprintf("  In network: %d | assigned=%d | grey=%d | %s\n",
              n_w_tot, n_w_asn, n_w_grey, kme_str))

  top_wrky <- gm_wrky[!is.na(gm_wrky$kME) & !is.na(gm_wrky$top_module) &
                      gm_wrky$top_module > 0L, ]
  top_wrky <- head(top_wrky[order(top_wrky$kME, decreasing = TRUE), ], 10L)
  if (nrow(top_wrky) > 0L) {
    cat("  Top-10 WRKY by kME:\n")
    for (i in seq_len(nrow(top_wrky))) {
      sym_w <- sym_lookup[top_wrky$gene_id[i]]
      cat(sprintf("    %s (%s)  module=%d  kME=%.4f\n",
                  top_wrky$gene_id[i],
                  ifelse(is.na(sym_w), "NA", sym_w),
                  top_wrky$top_module[i],
                  top_wrky$kME[i]))
    }
  }
} else if (length(wrky_ids) == 0L) {
  cat("  No WRKY IDs available.\n")
} else {
  cat("  Louvain failed — no WRKY readout.\n")
}

message(sprintf("\nStep 4 complete — feature plots + BON3/WRKY (%.1f min)",
                (proc.time() - t4)[["elapsed"]] / 60))

# ==============================================================================
# FINAL REPORT
# ==============================================================================

t_total <- (proc.time() - t_global)[["elapsed"]]

message("\n\n========== FINAL REPORT ==========")
message(sprintf("1. Network at |r|>=%.2f: n_pairs=%d  n_genes=%d  density=%.4f",
                MIN_ABS_R, n_pairs, n_genes, density))

if (!is.null(louvain_result)) {
  message(sprintf("2. Louvain: modules=%d | grey_rate=%.1f%% | size min/med/max=%s/%s/%s",
                  n_mods_l, 100 * n_grey_l / nrow(louvain_result$gene_module),
                  if (length(sizes_l) > 0L) min(sizes_l) else "NA",
                  if (length(sizes_l) > 0L) sprintf("%.0f", median(sizes_l)) else "NA",
                  if (length(sizes_l) > 0L) max(sizes_l) else "NA"))
} else {
  message("2. Louvain: FAILED")
}

if (!is.null(wgcna_result)) {
  message(sprintf("3. WGCNA: %s | modules=%d | grey_rate=%.1f%%",
                  wgcna_outcome, n_mods_w, 100 * n_grey_w / nrow(wgcna_result$gene_module)))
} else {
  message(sprintf("3. WGCNA: %s", wgcna_outcome))
}

message(sprintf("4. BON3 partners at |r|>=%.2f: %d", MIN_ABS_R, nrow(bon3_042)))
message(sprintf("5. WRKY IDs: %d", length(wrky_ids)))
message(sprintf("Total wall time: %.1f min", t_total / 60))
message("=== DONE ===")
