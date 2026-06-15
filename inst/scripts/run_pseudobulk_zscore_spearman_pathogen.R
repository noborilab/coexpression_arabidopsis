## FLAG-14 re-run: zscore_gene + Spearman + obs_subcluster(298 pts)
## Pareto-dominant design from Stage 1-2 observation-point sweep.
## Produces full analysis: per-condition networks, robustness, condition
## patterns, 4 module sets, module condition profiles, feature plots, BON3/WRKY.

suppressPackageStartupMessages({
  library(CoexprArabidopsis)
  library(igraph)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

t_global <- proc.time()

# ==============================================================================
# PARAMETERS
# ==============================================================================

SEURAT_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/",
  "Projects/SA_PTI_ETI_single_cell/",
  "SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed/",
  "combined_filtered.rds"
)

DATASET_ID  <- "pathogen_multiome"
CONDITIONS  <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")

# Condition-pattern labels: named vector mapping bit-pattern codes to human-readable
# labels for THIS dataset. Supply to characterize_condition_pattern(pattern_labels=).
# For a new dataset, replace with labels appropriate to your condition set.
# Patterns not listed here fall back to the generic "pattern_<bits>" label.
PATTERN_LABELS <- c(
  "0000" = "none",           "1111" = "constitutive_all",
  "1000" = "single_Mock",    "0100" = "single_DC3000",
  "0010" = "single_AvrRpt2", "0001" = "single_AvrRpm1",
  "0111" = "pan_pathogen",   "0011" = "ETI_shared"
)

STRATUM_VAR    <- "sample2"
STRATUM_LEVELS <- c(
  "00_Mock",
  "AvrRpm1_04h", "AvrRpm1_06h", "AvrRpm1_09h", "AvrRpm1_24h",
  "AvrRpt2_04h", "AvrRpt2_06h", "AvrRpt2_09h", "AvrRpt2_24h",
  "DC3000_04h",  "DC3000_06h",  "DC3000_09h",  "DC3000_24h"
)

SUBCLUSTER_COL <- "sub_clst_rna_20260610"  # dataset-specific: edit for new datasets

TF_META_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Nobori Lab (TSL) Team Folder/",
  "shared/datasets/from_Ben/for_tatsuya/data/motifs-2026/",
  "Athaliana_motifs_metadata.tsv"
)

RESULTS_DIR  <- file.path("results", DATASET_ID)
OUT_DIR      <- file.path(RESULTS_DIR, "pseudobulk_zscore_spearman")
OUT_MOD_DIR  <- file.path(OUT_DIR, "modules")
PLOTS_DIR    <- file.path(OUT_DIR, "plots")
SYMBOL_MAP_PATH <- file.path(RESULTS_DIR, "symbol_map.csv")
WRKY_PATH    <- file.path(RESULTS_DIR, "geneset_lookups", "WRKY_GGM_vs_PB.csv")

dir.create(OUT_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_MOD_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)

# Derive 4-level condition from sample2 column
.derive_condition <- function(vals) {
  cond <- sub("_(04h|06h|09h|24h).*", "", vals)
  cond[cond == "00_Mock"] <- "Mock"
  cond
}

# Module set configurations (same as run_official_modules_pathogen.R)
SET_CONFIGS <- list(
  large_wgcna  = list(method = "wgcna",   r_min = 0.5),
  large_louvain = list(method = "louvain", r_min = 0.5),
  small_wgcna  = list(method = "wgcna",   r_min = 0.6),
  small_louvain = list(method = "louvain", r_min = 0.6)
)

NAMED_LABELS <- c("constitutive_all", "pan_pathogen", "ETI_shared",
                  "single_Mock", "single_DC3000", "single_AvrRpt2",
                  "single_AvrRpm1", "none")

# ==============================================================================
# Load shared resources
# ==============================================================================

symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)
sym_lookup  <- setNames(symbol_map$gene_symbol, symbol_map$gene_id)
message("symbol_map: ", nrow(symbol_map), " entries")

.add_symbol <- function(df, col = "gene_id") {
  if (col %in% names(df)) df$gene_symbol <- sym_lookup[df[[col]]]
  df
}

# ==============================================================================
# Louvain module builder (same as run_official_modules_pathogen.R)
# ==============================================================================

.build_louvain_modules <- function(rob, r_score_min, min_module_size = 30L) {
  ps <- rob$pair_scores
  ps <- ps[!is.na(ps$R_score) & ps$R_score >= r_score_min, , drop = FALSE]
  if (nrow(ps) == 0L)
    stop("No pairs after r_score_min filter = ", r_score_min)

  edges <- data.frame(
    from   = ps$gene_id_A,
    to     = ps$gene_id_B,
    weight = pmin(abs(tanh(ps$z_bar)), 1.0),
    stringsAsFactors = FALSE
  )

  g  <- igraph::graph_from_data_frame(edges, directed = FALSE)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  memb      <- igraph::membership(cl)
  gene_ids  <- names(memb)
  top_lbl   <- as.integer(memb)

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
  if (is.null(hub_genes) || nrow(hub_genes) == 0L) {
    hub_genes <- data.frame(module_id = integer(), gene_id = character(),
                            gene_symbol = character(), kME = numeric(),
                            hub_rank = integer(), stringsAsFactors = FALSE)
  }

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

# Annotate and save one module set
.annotate_and_save_set <- function(set_name, mod_input, network_list, rob,
                                   tf_meta_path, out_dir) {
  t0      <- proc.time()
  set_dir <- file.path(out_dir, set_name)
  dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)

  message("  Annotating: ", set_name)

  mod_input$gene_module <- .add_symbol(mod_input$gene_module)
  mod_input$hub_genes   <- .add_symbol(mod_input$hub_genes)

  mod_input <- tryCatch(
    annotate_context(mod_input, network_list, ref_condition = "Mock"),
    error = function(e) {
      message("  annotate_context failed: ", conditionMessage(e)); mod_input
    }
  )

  mod_input <- tryCatch(
    annotate_go(mod_input, org_db = "org.At.tair.db", pval_cut = 0.05),
    error = function(e) {
      message("  annotate_go failed: ", conditionMessage(e)); mod_input
    }
  )

  if (file.exists(tf_meta_path)) {
    mod_input <- tryCatch(
      annotate_tfs(mod_input, tf_meta_path),
      error = function(e) {
        message("  annotate_tfs failed: ", conditionMessage(e)); mod_input
      }
    )
  } else {
    message("  TF file not found: ", tf_meta_path)
  }

  test_nets <- network_list[names(network_list) != "Mock"]
  pres <- tryCatch(
    compute_preservation_fallback(mod_input, test_nets),
    error = function(e) {
      message("  preservation failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(pres)) {
    mod_input$module_meta <- merge(mod_input$module_meta, pres,
                                   by = "module_id", all.x = TRUE)
    if ("zsummary.x" %in% names(mod_input$module_meta)) {
      mod_input$module_meta$zsummary           <- mod_input$module_meta$zsummary.y
      mod_input$module_meta$preservation_method <- mod_input$module_meta$preservation_method.y
      mod_input$module_meta$zsummary.x               <- NULL
      mod_input$module_meta$zsummary.y               <- NULL
      mod_input$module_meta$preservation_method.x    <- NULL
      mod_input$module_meta$preservation_method.y    <- NULL
    }
  }

  write.csv(mod_input$gene_module,  file.path(set_dir, "gene_module.csv"),  row.names = FALSE)
  write.csv(mod_input$module_meta,  file.path(set_dir, "module_meta.csv"),  row.names = FALSE)
  write.csv(mod_input$module_hier,  file.path(set_dir, "module_hier.csv"),  row.names = FALSE)
  write.csv(mod_input$hub_genes,    file.path(set_dir, "hub_genes.csv"),    row.names = FALSE)
  write.csv(mod_input$module_tfs,   file.path(set_dir, "module_tfs.csv"),   row.names = FALSE)
  write.csv(as.data.frame(mod_input$eigengenes),
            file.path(set_dir, "eigengenes.csv"), row.names = TRUE)
  saveRDS(mod_input, file.path(set_dir, "module_input.rds"))

  gm <- mod_input$gene_module
  message(sprintf(
    "  [%s] modules=%d | assigned=%d | grey=%d (%.1f%%) | TF_entries=%d | %.1f min",
    set_name,
    length(unique(gm$top_module[gm$top_module > 0L])),
    sum(gm$top_module > 0L, na.rm = TRUE),
    sum(is.na(gm$top_module) | gm$top_module == 0L),
    100 * sum(is.na(gm$top_module) | gm$top_module == 0L) / nrow(gm),
    nrow(mod_input$module_tfs),
    (proc.time() - t0)[["elapsed"]] / 60
  ))

  mod_input
}

# ==============================================================================
# STEP 1: Load Seurat once; build per-condition networks
# ==============================================================================

message("\n==== Step 1: Per-condition pseudobulk networks ====")
message("  Design: obs_subcluster(", SUBCLUSTER_COL, ") + zscore_gene + Spearman")
t1 <- proc.time()

bundle_full <- load_seurat(
  seurat_path    = SEURAT_PATH,
  dataset_id     = DATASET_ID,
  stratum_var    = STRATUM_VAR,
  stratum_levels = STRATUM_LEVELS,
  assay          = "RNA",
  slot           = "data",
  min_cells      = 10L,
  symbol_map     = symbol_map
)

# Derive 4-level condition; keep as character (not factor) to avoid subset issues
bundle_full$cell_meta$condition_4 <- .derive_condition(
  as.character(bundle_full$cell_meta[[STRATUM_VAR]])
)

message("Full bundle: ", nrow(bundle_full$counts), " genes x ",
        ncol(bundle_full$counts), " cells")
message("Condition cell counts:")
print(table(bundle_full$cell_meta$condition_4))

# Verify SUBCLUSTER_COL is present
if (!SUBCLUSTER_COL %in% names(bundle_full$cell_meta)) {
  stop("Column '", SUBCLUSTER_COL, "' not found in cell_meta. Available: ",
       paste(head(names(bundle_full$cell_meta), 20), collapse = ", "))
}

network_list <- list()

for (cond in CONDITIONS) {
  message("\n  -- Condition: ", cond, " --")

  idx <- bundle_full$cell_meta$condition_4 == cond
  if (!any(idx)) {
    warning("No cells for condition '", cond, "'; skipping.")
    next
  }

  sub_bundle <- list(
    counts     = bundle_full$counts[, idx, drop = FALSE],
    counts_raw = if (!is.null(bundle_full$counts_raw))
                   bundle_full$counts_raw[, idx, drop = FALSE] else NULL,
    cell_meta  = bundle_full$cell_meta[idx, , drop = FALSE],
    gene_meta  = bundle_full$gene_meta,
    stratum_spec = list(variable = "condition_4", levels = cond),
    dataset_id = bundle_full$dataset_id
  )

  obs <- obs_subcluster(sub_bundle, group_col = SUBCLUSTER_COL)
  obs$matrix <- normalize_obs(obs, method = "zscore_gene")

  message("  obs: ", ncol(obs$matrix), " subcluster points; computing Spearman...")
  t_cor <- proc.time()
  net <- coexpr_from_obs(obs, cor_type = "spearman", storage_cutoff = 0.1)
  message(sprintf("  Spearman done in %.1f min | %d edges stored",
                  (proc.time() - t_cor)[["elapsed"]] / 60,
                  nrow(net$edge_table)))

  # Save cor matrix
  cor_path <- file.path(OUT_DIR, paste0("cor_", cond, ".rds"))
  saveRDS(net$cor_mat, cor_path)
  message("  Saved: ", cor_path)

  # Assemble NetworkResult with n_pseudobulk for robustness layer
  network_list[[cond]] <- list(
    edge_table   = net$edge_table,
    gene_ids     = net$gene_ids,
    stratum_id   = cond,
    mode         = "pseudobulk_zscore_spearman",
    params       = list(n_pseudobulk = ncol(obs$matrix)),
    cor_type     = "spearman"
  )

  rm(obs, net, sub_bundle); gc()
}

rm(bundle_full); gc()
message(sprintf("Step 1 complete — %d per-condition networks built (%.1f min)",
                length(network_list),
                (proc.time() - t1)[["elapsed"]] / 60))

# ==============================================================================
# STEP 2: Robustness
# ==============================================================================

message("\n==== Step 2: Robustness ====")
t2 <- proc.time()

rob <- compute_robustness(network_list, k = 1.64, weight_cap = 30, fdr_method = "BH")

n_full <- nrow(rob$pair_scores)
n_r05  <- sum(rob$pair_scores$R_score >= 0.5, na.rm = TRUE)
n_r06  <- sum(rob$pair_scores$R_score >= 0.6, na.rm = TRUE)
message(sprintf("  Pairs: full=%d | R>=0.5: %d | R>=0.6: %d", n_full, n_r05, n_r06))

write.csv(rob$pair_scores,
          file.path(OUT_DIR, "pair_scores_full.csv"), row.names = FALSE)
write.csv(rob$pair_scores[!is.na(rob$pair_scores$R_score) & rob$pair_scores$R_score >= 0.5, ],
          file.path(OUT_DIR, "pair_scores_r0.5.csv"), row.names = FALSE)
write.csv(rob$pair_scores[!is.na(rob$pair_scores$R_score) & rob$pair_scores$R_score >= 0.6, ],
          file.path(OUT_DIR, "pair_scores_r0.6.csv"), row.names = FALSE)
saveRDS(rob, file.path(OUT_DIR, "robustness_result.rds"))

message(sprintf("Step 2 complete — pair_scores_full.csv: %d pairs (%.1f min)",
                n_full, (proc.time() - t2)[["elapsed"]] / 60))

# ==============================================================================
# STEP 3: Condition-pattern characterisation
# ==============================================================================

message("\n==== Step 3: Condition-pattern characterisation ====")
t3 <- proc.time()

cp <- characterize_condition_pattern(rob, network_list, condition_order = CONDITIONS,
                                     pattern_labels = PATTERN_LABELS)
write.csv(cp, file.path(OUT_DIR, "pair_condition_patterns.csv"), row.names = FALSE)

pat_counts <- sort(table(cp$pattern_label), decreasing = TRUE)
message("  Pattern distribution:")
for (nm in names(pat_counts))
  message(sprintf("    %-24s %d", nm, pat_counts[[nm]]))

message(sprintf("Step 3 complete — pair_condition_patterns.csv: %d pairs (%.1f min)",
                nrow(cp), (proc.time() - t3)[["elapsed"]] / 60))

# ==============================================================================
# STEP 4: Module construction (4 sets)
# ==============================================================================

message("\n==== Step 4: Module construction ====")
t4 <- proc.time()

mod_results <- list()

for (set_name in names(SET_CONFIGS)) {
  cfg <- SET_CONFIGS[[set_name]]
  message(sprintf("\n  === SET: %s (method=%s, R>=%s) ===",
                  set_name, cfg$method, cfg$r_min))
  t_set <- proc.time()

  tryCatch({
    mod_input <- if (cfg$method == "wgcna") {
      build_wgcna_modules(
        rob             = rob,
        network_list    = network_list,
        r_score_min     = cfg$r_min,
        soft_power      = 1L,
        merge_cut       = 0.25,
        min_module_size = 30L,
        sub_merge_cut   = 0.10
      )
    } else {
      .build_louvain_modules(rob, r_score_min = cfg$r_min, min_module_size = 30L)
    }

    mod_input$method            <- if (cfg$method == "wgcna") "wgcna_p1" else "louvain"
    mod_input$graph             <- if (cfg$r_min == 0.5) "large" else "small"
    mod_input$r_score_threshold <- cfg$r_min

    gm        <- mod_input$gene_module
    n_modules <- length(unique(gm$top_module[gm$top_module > 0L]))
    n_grey    <- sum(is.na(gm$top_module) | gm$top_module == 0L)
    message(sprintf("  Build done: %d modules | %d grey (%.1f%%) | %.1f min",
                    n_modules, n_grey,
                    100 * n_grey / nrow(gm),
                    (proc.time() - t_set)[["elapsed"]] / 60))

    mod_input <- .annotate_and_save_set(
      set_name      = set_name,
      mod_input     = mod_input,
      network_list  = network_list,
      rob           = rob,
      tf_meta_path  = TF_META_PATH,
      out_dir       = OUT_MOD_DIR
    )

    mod_results[[set_name]] <- mod_input

  }, error = function(e) {
    message("  ERROR in set '", set_name, "': ", conditionMessage(e))
  })
}

message(sprintf("Step 4 complete — %d module sets built (%.1f min)",
                length(mod_results), (proc.time() - t4)[["elapsed"]] / 60))

# ==============================================================================
# STEP 5: Module condition-pattern profiles
# ==============================================================================

message("\n==== Step 5: Module condition-pattern profiles ====")
t5 <- proc.time()

all_labels  <- sort(unique(cp$pattern_label))
mixed_lbls  <- sort(grep("^mixed_", all_labels, value = TRUE))
label_order <- c(NAMED_LABELS[NAMED_LABELS %in% all_labels], mixed_lbls)

all_mod_rows <- list()

for (set_name in names(SET_CONFIGS)) {
  gm_path <- file.path(OUT_MOD_DIR, set_name, "gene_module.csv")
  if (!file.exists(gm_path)) {
    message("  ", set_name, ": gene_module.csv not found; skipping.")
    next
  }
  gm <- read.csv(gm_path, stringsAsFactors = FALSE)
  gm <- gm[gm$top_module != 0L, , drop = FALSE]

  gene_to_mod <- setNames(gm$top_module, gm$gene_id)
  mod_A <- gene_to_mod[cp$gene_id_A]
  mod_B <- gene_to_mod[cp$gene_id_B]
  keep  <- !is.na(mod_A) & !is.na(mod_B) & (mod_A == mod_B)
  intra <- cp[keep, , drop = FALSE]
  intra$module_id <- as.integer(mod_A[keep])

  message(sprintf("  %s: %d genes in %d modules | %d intra-module pairs",
                  set_name, nrow(gm), length(unique(gm$top_module)),
                  nrow(intra)))

  mods <- sort(unique(intra$module_id))
  rows <- lapply(mods, function(mid) {
    sub <- intra[intra$module_id == mid, ]
    n   <- nrow(sub)
    pat_tbl  <- table(sub$pattern_label)
    pat_frac <- pat_tbl / n
    dom_pat  <- names(which.max(pat_frac))
    dom_frac <- as.numeric(max(pat_frac))

    frac_cols <- setNames(
      vapply(label_order, function(lbl)
        if (lbl %in% names(pat_frac)) as.numeric(pat_frac[[lbl]]) else 0.0,
        numeric(1)),
      paste0("frac_", label_order)
    )

    c(list(
      set = set_name, module_id = mid, n_intra_edges = n,
      dominant_pattern = dom_pat, dominant_fraction = dom_frac,
      w_Mock = mean(sub$w_Mock, na.rm = TRUE),
      w_DC3000 = mean(sub$w_DC3000, na.rm = TRUE),
      w_AvrRpt2 = mean(sub$w_AvrRpt2, na.rm = TRUE),
      w_AvrRpm1 = mean(sub$w_AvrRpm1, na.rm = TRUE),
      module_specificity_index = mean(sub$specificity_index, na.rm = TRUE),
      n_conditions_active_mean = mean(sub$n_conditions_active, na.rm = TRUE)
    ), as.list(frac_cols))
  })

  set_df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))

  per_set_path <- file.path(OUT_MOD_DIR, set_name, "module_condition_patterns.csv")
  write.csv(set_df, per_set_path, row.names = FALSE)
  message("  Saved: ", per_set_path)

  all_mod_rows[[set_name]] <- set_df
}

combined_df   <- do.call(rbind, all_mod_rows)
combined_path <- file.path(OUT_MOD_DIR, "all_modules_condition_patterns.csv")
write.csv(combined_df, combined_path, row.names = FALSE)
message("Saved: ", combined_path, " (", nrow(combined_df), " module rows)")

message(sprintf("Step 5 complete — module condition profiles for %d sets (%.1f min)",
                length(all_mod_rows), (proc.time() - t5)[["elapsed"]] / 60))

# ==============================================================================
# STEP 6: Feature plots
# ==============================================================================

message("\n==== Step 6: Feature plots ====")
t6 <- proc.time()

suppressPackageStartupMessages(library(Seurat))

PURPLE_COLS <- c("lightgray", "#BFD3E6", "#9EBCDA", "#8C96C6",
                  "#8C6BB1", "#88419D", "#810F7C", "#4D004B")
TOP_N <- 4L

message("  Loading Seurat object for feature plots...")
sobj <- readRDS(SEURAT_PATH)

# Derive condition column (matches featureplot_modules.R)
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
    { sym2 <- symbol_map$gene_symbol[symbol_map$gene_id == gene_id]
      if (length(sym2) > 0 && !is.na(sym2[1])) sym2[1] else character(0) }
  )
  for (cand in cands) if (cand %in% RNA_GENES) return(cand)
  NA_character_
}

.make_label <- function(gene_id, gene_symbol) {
  if (!is.na(gene_symbol) && nchar(trimws(gene_symbol)) > 0)
    paste0(gene_id, " / ", gene_symbol)
  else
    gene_id
}

for (set_name in names(SET_CONFIGS)) {
  hub_path <- file.path(OUT_MOD_DIR, set_name, "hub_genes.csv")
  gm_path  <- file.path(OUT_MOD_DIR, set_name, "gene_module.csv")
  if (!file.exists(hub_path) || !file.exists(gm_path)) {
    message("  ", set_name, ": missing files; skipping feature plots.")
    next
  }

  hub_genes_df <- read.csv(hub_path, stringsAsFactors = FALSE)
  gene_mod_df  <- read.csv(gm_path,  stringsAsFactors = FALSE)

  all_modules <- sort(unique(hub_genes_df$module_id))
  pdf_path    <- file.path(PLOTS_DIR, paste0("featureplots_", set_name, ".pdf"))

  message(sprintf("  %s: %d modules → %s", set_name, length(all_modules), pdf_path))

  pdf(pdf_path,
      width  = length(CONDITIONS) * 5,
      height = (TOP_N + 1L) * 5,
      onefile = TRUE)

  for (mid in all_modules) {
    top_genes <- hub_genes_df[hub_genes_df$module_id == mid, ]
    top_genes <- top_genes[order(top_genes$hub_rank), ]
    top_genes <- head(top_genes, TOP_N)

    if (nrow(top_genes) == 0L) next

    top_genes$feature <- mapply(
      .resolve_feature, top_genes$gene_id,
      ifelse(is.na(top_genes$gene_symbol), NA_character_, top_genes$gene_symbol),
      SIMPLIFY = TRUE
    )
    top_genes_ok <- top_genes[!is.na(top_genes$feature), , drop = FALSE]
    if (nrow(top_genes_ok) == 0L) next

    # Module score
    all_ids  <- gene_mod_df$gene_id[gene_mod_df$top_module == mid]
    all_rna  <- vapply(all_ids, function(gid) {
      sym2 <- symbol_map$gene_symbol[symbol_map$gene_id == gid]
      .resolve_feature(gid, if (length(sym2) > 0 && !is.na(sym2[1])) sym2[1] else NA_character_)
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

        fp <- FeaturePlot(sobj,
                          features   = feat,
                          split.by   = "condition",
                          reduction  = "umap",
                          pt.size    = 0.5,
                          order      = TRUE,
                          max.cutoff = "q99",
                          min.cutoff = "q1")
        fp <- fp &
          scale_colour_gradientn(colours = PURPLE_COLS, na.value = "lightgray") &
          NoAxes()
        for (k in seq_along(fp)) {
          new_title <- if (k == 1) paste0(label, "\n(", CONDITIONS[k], ")") else CONDITIONS[k]
          fp[[k]] <- fp[[k]] + labs(title = new_title) +
            theme(plot.title = element_text(size = 9, face = "bold"))
        }
        fp
      })

      score_panels <- lapply(seq_along(CONDITIONS), function(k) {
        cond_k <- CONDITIONS[k]
        df_k   <- umap_df[umap_df$condition == cond_k, ]
        df_k   <- df_k[order(df_k$module_score, na.last = FALSE), ]
        title  <- if (k == 1) paste0("M", mid, " score\n(", cond_k, ")") else cond_k

        ggplot(df_k, aes(x = UMAP_1, y = UMAP_2, color = module_score)) +
          geom_point(size = 0.5) +
          scale_color_gradientn(colours = PURPLE_COLS, na.value = "lightgray",
                                name = "score") +
          labs(title = title) +
          NoAxes() +
          theme(plot.background  = element_blank(),
                panel.background = element_blank(),
                plot.title       = element_text(size = 9, face = "bold"))
      })
      score_row <- wrap_plots(score_panels, nrow = 1)
      final_fig <- wrap_plots(c(gene_rows, list(score_row)), ncol = 1)
      print(final_fig)

    }, error = function(e) {
      message(sprintf("  ERROR M%s (%s): %s", mid, set_name, conditionMessage(e)))
    })
  }

  dev.off()
  message(sprintf("  Saved: %s", pdf_path))
}

rm(sobj); gc()
message(sprintf("Step 6 complete — feature plots (%.1f min)",
                (proc.time() - t6)[["elapsed"]] / 60))

# ==============================================================================
# STEP 7: BON3 / WRKY post-hoc sanity readout
# ==============================================================================

message("\n==== Step 7: BON3 / WRKY sanity readout ====")

BON3_ID <- "AT1G08860"

# Identify WRKY AT-IDs
wrky_ids <- character(0)
if (file.exists(WRKY_PATH)) {
  wrky_df  <- read.csv(WRKY_PATH, stringsAsFactors = FALSE)
  wrky_ids <- wrky_df$gene_id
  message("  WRKY IDs loaded from WRKY_GGM_vs_PB.csv: ", length(wrky_ids))
} else if (file.exists(TF_META_PATH)) {
  tf_meta  <- read.table(TF_META_PATH, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, quote = "")
  fam_col  <- if ("tf_family" %in% names(tf_meta)) "tf_family"
               else if ("family" %in% names(tf_meta)) "family"
               else if ("class" %in% names(tf_meta)) "class"
               else NA_character_
  if (!is.na(fam_col)) {
    wrky_rows <- tf_meta[grepl("WRKY", tf_meta[[fam_col]], ignore.case = TRUE), ]
    wrky_ids  <- wrky_rows$motif_id
    message("  WRKY IDs from TF metadata: ", length(wrky_ids))
  }
} else {
  message("  WARNING: WRKY source not found; WRKY section will be empty.")
}

# BON3 partners from robustness
ps <- rob$pair_scores
bon3_ps <- ps[ps$gene_id_A == BON3_ID | ps$gene_id_B == BON3_ID, , drop = FALSE]
if (nrow(bon3_ps) > 0L) {
  bon3_ps$partner <- ifelse(bon3_ps$gene_id_A == BON3_ID,
                            bon3_ps$gene_id_B, bon3_ps$gene_id_A)
  bon3_ps <- bon3_ps[order(bon3_ps$R_score, decreasing = TRUE), ]
}

cat("\n")
cat("=== BON3 / WRKY Summary ===\n\n")
cat(sprintf("BON3 (%s) across 4 module sets:\n", BON3_ID))
cat(sprintf("  Partners in pair_scores (any R_score): %d\n", nrow(bon3_ps)))
cat(sprintf("  Partners at R_score >= 0.5: %d\n",
            sum(bon3_ps$R_score >= 0.5, na.rm = TRUE)))

header_fmt <- "  %-20s %12s %10s %14s\n"
row_fmt    <- "  %-20s %12s %10s %14s\n"
cat(sprintf(header_fmt, "Module set", "module_id", "kME", "n_partners(R>=0)"))
cat(sprintf(header_fmt, "----------", "---------", "---", "---------------"))

for (set_name in names(SET_CONFIGS)) {
  mi <- mod_results[[set_name]]
  if (is.null(mi)) { cat(sprintf(row_fmt, set_name, "FAILED", "-", "-")); next }
  gm <- mi$gene_module
  r  <- gm[gm$gene_id == BON3_ID, , drop = FALSE]
  if (nrow(r) == 0L) {
    cat(sprintf(row_fmt, set_name, "not in network", "NA", "-"))
    next
  }
  top_mod <- r$top_module[1L]
  kme_val <- r$kME[1L]
  mod_lbl <- if (!is.na(top_mod) && top_mod == 0L) "grey(0)" else as.character(top_mod)
  cat(sprintf(row_fmt, set_name, mod_lbl,
              sprintf("%.5f", kme_val), nrow(bon3_ps)))
}

if (nrow(bon3_ps) > 0L && length(sym_lookup) > 0L) {
  bon3_top <- head(bon3_ps, 10L)
  bon3_top$partner_symbol <- sym_lookup[bon3_top$partner]
  cat("\nBON3 top-10 partners (by R_score):\n")
  for (i in seq_len(nrow(bon3_top))) {
    sym_str <- ifelse(is.na(bon3_top$partner_symbol[i]), "NA",
                      bon3_top$partner_symbol[i])
    cat(sprintf("  %2d. %s (%s)  R_score=%.4f  z_bar=%.4f\n",
                i, bon3_top$partner[i], sym_str,
                bon3_top$R_score[i], bon3_top$z_bar[i]))
  }
}

cat(sprintf("\nWRKY family (%d AT-IDs) across 4 module sets:\n", length(wrky_ids)))

if (length(wrky_ids) > 0L) {
  for (set_name in names(SET_CONFIGS)) {
    mi <- mod_results[[set_name]]
    if (is.null(mi)) { cat(sprintf("  %s: FAILED\n", set_name)); next }
    gm <- mi$gene_module
    gm_wrky <- gm[gm$gene_id %in% wrky_ids, , drop = FALSE]

    if (nrow(gm_wrky) == 0L) {
      cat(sprintf("  %s: 0 WRKY genes in network\n", set_name))
      next
    }

    n_total    <- nrow(gm_wrky)
    n_assigned <- sum(gm_wrky$top_module > 0L, na.rm = TRUE)
    n_grey     <- n_total - n_assigned

    kme_assigned <- gm_wrky$kME[!is.na(gm_wrky$kME) & gm_wrky$top_module > 0L]
    kme_stats <- if (length(kme_assigned) > 0L) {
      sprintf("kME: min=%.3f med=%.3f max=%.3f",
              min(kme_assigned), median(kme_assigned), max(kme_assigned))
    } else "kME: N/A"

    cat(sprintf("  %s: %d in network | assigned=%d | grey=%d | %s\n",
                set_name, n_total, n_assigned, n_grey, kme_stats))

    # Top-5 WRKY by kME
    top_wrky <- gm_wrky[!is.na(gm_wrky$kME) & gm_wrky$top_module > 0L, ]
    top_wrky <- head(top_wrky[order(top_wrky$kME, decreasing = TRUE), ], 5L)
    if (nrow(top_wrky) > 0L) {
      cat("    Top-5 by kME:\n")
      for (i in seq_len(nrow(top_wrky))) {
        sym_w <- sym_lookup[top_wrky$gene_id[i]]
        cat(sprintf("      %s (%s) mod=%d kME=%.4f\n",
                    top_wrky$gene_id[i],
                    ifelse(is.na(sym_w), "NA", sym_w),
                    top_wrky$top_module[i],
                    top_wrky$kME[i]))
      }
    }
  }
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

t_total <- (proc.time() - t_global)[["elapsed"]]

message("\n\n========== FINAL REPORT ==========")
message(sprintf("1. pair_scores_full.csv: %d pairs", n_full))
message(sprintf("2. R_score >= 0.5: %d pairs | R_score >= 0.6: %d pairs", n_r05, n_r06))
message("3. Module sets:")
for (set_name in names(SET_CONFIGS)) {
  mi <- mod_results[[set_name]]
  if (is.null(mi)) {
    message(sprintf("   %-20s : FAILED", set_name)); next
  }
  gm     <- mi$gene_module
  n_mods <- length(unique(gm$top_module[gm$top_module > 0L]))
  n_grey <- sum(is.na(gm$top_module) | gm$top_module == 0L)
  n_tot  <- nrow(gm)
  message(sprintf("   %-20s : %d modules | grey=%d/%d (%.1f%%)",
                  set_name, n_mods, n_grey, n_tot, 100 * n_grey / n_tot))
}
message(sprintf("4. BON3 (%s): %d partners (any R_score)", BON3_ID, nrow(bon3_ps)))
message(sprintf("5. WRKY IDs available: %d", length(wrky_ids)))
message(sprintf("Total wall time: %.1f min", t_total / 60))
message("=== DONE ===")
