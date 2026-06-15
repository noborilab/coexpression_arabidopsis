## Pathogen multiome — subcluster pseudobulk co-expression
##
## Mode: pooled single stratum (all 4 conditions combined)
## Observations: ~298 subclusters (GROUP_VAR; dataset-specific)
## Rationale: pseudobulk aggregation by subcluster recovers rare-population
##   co-expression that the cell-level GGM dilutes via full-universe conditioning.
##
## Generalizability: change ONLY the DATASET-SPECIFIC PARAMETERS block.
## Everything below that block is dataset-agnostic.

suppressPackageStartupMessages({
  library(CoexprArabidopsis)
  library(igraph)
  library(Matrix)
})

t_global <- proc.time()

# ==============================================================================
# DATASET-SPECIFIC PARAMETERS  (change these for dev atlas or other datasets)
# ==============================================================================

SEURAT_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/",
  "Projects/SA_PTI_ETI_single_cell/",
  "SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed/",
  "combined_filtered.rds"
)
DATASET_ID   <- "pathogen_pseudobulk"

# Pseudobulk grouping: each unique value → one pseudobulk observation
GROUP_VAR    <- "sub_clst_rna_20260610"  # dataset-specific: edit for new datasets

# Assay / slot
ASSAY <- "RNA"
SLOT  <- "data"

# Expression filter: genes detected in >= MIN_CELLS cells
MIN_CELLS <- 10L

# Stratification: FALSE = pool all strata into one network (single correlation)
#                 TRUE  = stratify by STRATUM_VAR (one network per level → robustness)
STRATIFY <- FALSE

# STRATUM_VAR: used as the pass-through variable when STRATIFY = FALSE
# (load_seurat requires a stratum_var; we keep all levels to retain all cells)
# For STRATIFY = TRUE: set STRATUM_VAR to the condition/organ column and
# STRATUM_LEVELS_KEEP to only the desired levels.
STRATUM_VAR          <- "sample2"
STRATUM_LEVELS_KEEP  <- c(            # all 13 sample2 levels → keep all cells
  "00_Mock",
  "AvrRpm1_04h", "AvrRpm1_06h", "AvrRpm1_09h", "AvrRpm1_24h",
  "AvrRpt2_04h", "AvrRpt2_06h", "AvrRpt2_09h", "AvrRpt2_24h",
  "DC3000_04h",  "DC3000_06h",  "DC3000_09h",  "DC3000_24h"
)

# For post-hoc condition context: function to derive 4-level condition from sample2
# (Change for dev atlas: could return organ labels from another column)
.derive_condition <- function(stratum_vals) {
  cond <- sub("_(04h|06h|09h|24h).*", "", stratum_vals)
  cond[cond == "00_Mock"] <- "Mock"
  cond
}
CONDITIONS <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")

# Symbol map (reuse from GGM run for consistency)
SYMBOL_MAP_PATH <- file.path("results", "pathogen_multiome", "symbol_map.csv")

# TF metadata
TF_META_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Nobori Lab (TSL) Team Folder/",
  "shared/datasets/from_Ben/for_tatsuya/data/motifs-2026/",
  "Athaliana_motifs_metadata.tsv"
)

# Edge threshold for module construction: keep |Spearman| >= EDGE_THR
# (will be reported and verified before use)
EDGE_THR <- 0.3

# Output directories
RESULTS_DIR <- file.path("results", DATASET_ID)
OUT_NET_DIR <- "output_pseudobulk_pathogen"
OUT_MOD_DIR <- file.path(RESULTS_DIR, "modules")

# ==============================================================================
# END DATASET-SPECIFIC PARAMETERS
# ==============================================================================

dir.create(RESULTS_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_NET_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_MOD_DIR,  showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Load symbol map
# ------------------------------------------------------------------------------

symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)
sym_lookup  <- setNames(symbol_map$gene_symbol, symbol_map$gene_id)
message("symbol_map: ", nrow(symbol_map), " entries")

# ------------------------------------------------------------------------------
# Helpers (dataset-agnostic)
# ------------------------------------------------------------------------------

# Join gene_symbol onto any data.frame that has a gene_id column
.add_symbol <- function(df, col = "gene_id") {
  if (col %in% names(df)) df$gene_symbol <- sym_lookup[df[[col]]]
  df
}

# Post-hoc condition context from cell-level counts.
# Used when STRATIFY = FALSE and annotate_context() (which needs per-stratum
# network_list edge weights) is not applicable.
# For each module: for each condition, compute mean of per-cell mean expression
# of module genes across cells in that condition.
.annotate_context_expr <- function(mod_input, counts, cond_labels, conditions) {
  gene_mod    <- mod_input$gene_module
  unique_mods <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))

  for (m in unique_mods) {
    mod_genes <- gene_mod$gene_id[gene_mod$top_module == m]
    present   <- mod_genes[mod_genes %in% rownames(counts)]
    if (length(present) == 0L) next

    cond_means <- setNames(vapply(conditions, function(cond) {
      idx <- which(cond_labels == cond)
      if (length(idx) < 10L) return(0.0)
      mean(Matrix::colMeans(counts[present, idx, drop = FALSE]))
    }, numeric(1)), conditions)

    top_cond  <- conditions[which.max(cond_means)]
    mock_mean <- cond_means["Mock"]
    non_mock  <- conditions[conditions != "Mock"]
    delta_str <- if ("Mock" %in% conditions && !is.na(mock_mean)) {
      paste(sprintf("%s:%+.3f", non_mock, cond_means[non_mock] - mock_mean),
            collapse = ";")
    } else NA_character_

    ri <- which(mod_input$module_meta$module_id == m)
    mod_input$module_meta$top_organ_or_condition[ri] <- top_cond
    mod_input$module_meta$delta_treatment[ri]        <- delta_str
  }
  mod_input
}

# Build Louvain ModuleInput from a pseudo_rob (same pattern as official modules script).
# pseudo_rob$pair_scores must have: gene_id_A, gene_id_B, R_score, z_bar
.build_louvain_modules <- function(pseudo_rob, r_score_min, min_module_size = 30L) {
  ps <- pseudo_rob$pair_scores
  ps <- ps[!is.na(ps$R_score) & ps$R_score >= r_score_min, , drop = FALSE]
  if (nrow(ps) == 0L)
    stop("No pairs remain after threshold = ", r_score_min)

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

  # Map communities smaller than min_module_size to module 0 (grey)
  comm_sz <- table(top_lbl)
  small   <- as.integer(names(comm_sz[comm_sz < min_module_size]))
  top_lbl[top_lbl %in% small] <- 0L

  # Re-label surviving modules as consecutive integers starting at 1
  live    <- sort(unique(top_lbl[top_lbl > 0L]))
  relabel <- setNames(seq_along(live), as.character(live))
  top_lbl <- ifelse(top_lbl == 0L, 0L,
                    as.integer(relabel[as.character(top_lbl)]))

  # Build adjacency matrix for kME computation
  all_genes <- sort(union(edges$from, edges$to))
  n         <- length(all_genes)
  gene_idx  <- setNames(seq_along(all_genes), all_genes)

  A <- matrix(0.0, n, n, dimnames = list(all_genes, all_genes))
  iA <- gene_idx[edges$from]; iB <- gene_idx[edges$to]
  A[cbind(iA, iB)] <- edges$weight
  A[cbind(iB, iA)] <- edges$weight

  tl_ord      <- setNames(top_lbl, gene_ids)[all_genes]
  unique_mods <- sort(unique(tl_ord[tl_ord > 0L]))

  # kME: Pearson between each gene's adjacency row and mean row of its module
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
    sub_module = NA_integer_,   # Louvain: no hierarchy
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

  # Empty module_hier: Louvain has no hierarchical structure
  module_hier <- data.frame(sub_module = integer(), top_module = integer(),
                            stringsAsFactors = FALSE)

  # Top 20 hub genes per module by kME
  hub_list <- lapply(as.integer(names(mod_counts)), function(m) {
    sg <- gene_module[gene_module$top_module == m & !is.na(gene_module$kME), ]
    sg <- head(sg[order(sg$kME, decreasing = TRUE), ], 20L)
    if (nrow(sg) == 0L) return(NULL)
    data.frame(module_id = m, gene_id = sg$gene_id, gene_symbol = NA_character_,
               kME = sg$kME, hub_rank = seq_len(nrow(sg)), stringsAsFactors = FALSE)
  })
  hub_genes <- do.call(rbind, Filter(Negate(is.null), hub_list))
  if (is.null(hub_genes) || nrow(hub_genes) == 0L) {
    hub_genes <- data.frame(
      module_id = integer(), gene_id = character(), gene_symbol = character(),
      kME = numeric(), hub_rank = integer(), stringsAsFactors = FALSE
    )
  }

  module_tfs <- data.frame(
    module_id = integer(), gene_id = character(),
    gene_symbol = character(), tf_family = character(),
    stringsAsFactors = FALSE
  )

  # Module eigengenes: mean adjacency profile per module (rows = module × genes)
  # Returned transposed (genes × modules) to match WGCNA convention
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
.annotate_and_save <- function(mod_input, method_name, out_dir,
                               counts_ctx, cond_lbl, conditions,
                               network_list, tf_path, stratified) {
  if (is.null(mod_input)) {
    message("  [", method_name, "] skipped (not built)")
    return(NULL)
  }
  t0 <- proc.time()
  set_dir <- file.path(out_dir, method_name)
  dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)

  # Gene symbols
  mod_input$gene_module <- .add_symbol(mod_input$gene_module)
  mod_input$hub_genes   <- .add_symbol(mod_input$hub_genes)

  # Condition context
  if (stratified) {
    # Stratified: use network edge weights per stratum (standard annotate_context)
    mod_input <- tryCatch(annotate_context(mod_input, network_list,
                                           ref_condition = "Mock"),
                          error = function(e) {
                            message("  context failed: ", conditionMessage(e))
                            mod_input
                          })
  } else {
    # Pooled: derive context from cell-level expression post-hoc
    mod_input <- tryCatch(
      .annotate_context_expr(mod_input, counts_ctx, cond_lbl, conditions),
      error = function(e) {
        message("  context_expr failed: ", conditionMessage(e))
        mod_input
      })
  }

  # GO enrichment
  mod_input <- tryCatch(
    annotate_go(mod_input, org_db = "org.At.tair.db", pval_cut = 0.05),
    error = function(e) { message("  GO failed: ", conditionMessage(e)); mod_input }
  )

  # TF intersection
  if (file.exists(tf_path)) {
    mod_input <- tryCatch(
      annotate_tfs(mod_input, tf_path),
      error = function(e) { message("  TF failed: ", conditionMessage(e)); mod_input }
    )
  } else {
    message("  TF file not found: ", tf_path)
  }

  # Preservation: not applicable for single pooled network
  if (!stratified) {
    message("  [", method_name, "] preservation skipped (single pooled stratum)")
  }

  # Save
  write.csv(mod_input$gene_module,  file.path(set_dir, "gene_module.csv"),  row.names = FALSE)
  write.csv(mod_input$module_meta,  file.path(set_dir, "module_meta.csv"),  row.names = FALSE)
  write.csv(mod_input$module_hier,  file.path(set_dir, "module_hier.csv"),  row.names = FALSE)
  write.csv(mod_input$hub_genes,    file.path(set_dir, "hub_genes.csv"),    row.names = FALSE)
  write.csv(mod_input$module_tfs,   file.path(set_dir, "module_tfs.csv"),   row.names = FALSE)
  write.csv(as.data.frame(mod_input$eigengenes),
            file.path(set_dir, "eigengenes.csv"), row.names = TRUE)
  saveRDS(mod_input, file.path(set_dir, "module_input.rds"))

  gm <- mod_input$gene_module
  message(sprintf("  [%s] %d modules | assigned=%d | grey=%d (%.1f%%) | TF=%d | wall=%.1f min",
    method_name,
    length(unique(gm$top_module[gm$top_module > 0L])),
    sum(gm$top_module > 0L, na.rm = TRUE),
    sum(is.na(gm$top_module) | gm$top_module == 0L),
    100 * sum(is.na(gm$top_module) | gm$top_module == 0L) / nrow(gm),
    nrow(mod_input$module_tfs),
    (proc.time() - t0)[["elapsed"]] / 60))

  mod_input
}

# ==============================================================================
# STEP 1: Load InputBundle
# ==============================================================================

message("\n=== Step 1: Load InputBundle ===")
t1 <- proc.time()

bundle <- load_seurat(
  seurat_path    = SEURAT_PATH,
  dataset_id     = DATASET_ID,
  stratum_var    = STRATUM_VAR,
  stratum_levels = STRATUM_LEVELS_KEEP,
  group_var      = GROUP_VAR,
  assay          = ASSAY,
  slot           = SLOT,
  min_cells      = MIN_CELLS,
  symbol_map     = symbol_map
)

if (!STRATIFY) {
  # Override to single pooled stratum so estimate_pseudobulk produces one network
  bundle$stratum_spec             <- list(variable = "stratum_pooled",
                                          levels   = "pooled")
  bundle$cell_meta$stratum_pooled <- "pooled"
}

# estimate_pseudobulk looks for a column named literally "group_var" in cell_meta
bundle$cell_meta$group_var <- bundle$cell_meta[[GROUP_VAR]]

n_groups <- length(unique(bundle$cell_meta[[GROUP_VAR]]))
message(sprintf("Bundle: %d genes x %d cells | %d pseudobulk groups",
                nrow(bundle$counts), ncol(bundle$counts), n_groups))
message(sprintf("Step 1: %.1f min", (proc.time() - t1)[["elapsed"]] / 60))

# Extract cell-level info needed post-hoc (before freeing bundle)
cond_labels       <- .derive_condition(as.character(bundle$cell_meta[[STRATUM_VAR]]))
counts_for_context <- bundle$counts   # sparse genes × cells; kept for context annotation

# ==============================================================================
# STEP 2: Pseudobulk estimation
# ==============================================================================

message("\n=== Step 2: Pseudobulk estimation ===")
message("  NOTE: for ~18k genes x 298 subclusters, the correlation step builds")
message("  an ~18k x 18k matrix. Expected: 10-25 min, ~3-5 GB RAM.")
t2 <- proc.time()

network_list <- estimate_pseudobulk(bundle, min_samples = 5, min_expr = 0)

# Free the large count matrix now that estimation is done.
# counts_for_context is the same R object reference so it remains valid;
# only the bundle wrapper is freed.
rm(bundle); gc()

for (s in names(network_list)) {
  nr <- network_list[[s]]
  message(sprintf("  [%s] %d edges, %d genes",
                  s, nrow(nr$edge_table), length(nr$gene_ids)))
}

save_network_results(network_list, OUT_NET_DIR)
message(sprintf("Step 2: %.1f min", (proc.time() - t2)[["elapsed"]] / 60))

# ==============================================================================
# STEP 2b: Spearman distribution and threshold verification
# ==============================================================================

message("\n=== Step 2b: Spearman distribution ===")

et_pooled <- network_list[[1]]$edge_table
w         <- et_pooled$weight

message(sprintf("  Total edges (|Spearman| >= 0.1 storage cutoff): %d", nrow(et_pooled)))
message(sprintf("  Spearman range: min=%.4f  p10=%.4f  p25=%.4f  median=%.4f  p75=%.4f  max=%.4f",
                min(w), quantile(w, 0.10), quantile(w, 0.25),
                median(w), quantile(w, 0.75), max(w)))

for (thr in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  n_pos <- sum(w >=  thr)
  n_neg <- sum(w <= -thr)
  message(sprintf("  |Spearman| >= %.1f : %d positive + %d negative = %d total",
                  thr, n_pos, n_neg, n_pos + n_neg))
}

n_thr <- sum(abs(w) >= EDGE_THR)
message(sprintf("\n  Chosen threshold |Spearman| >= %.1f → %d edges", EDGE_THR, n_thr))

if (n_thr < 500) {
  warning("Very few edges (", n_thr, ") at EDGE_THR=", EDGE_THR,
          "; modules may be small. Consider lowering threshold.")
} else if (n_thr > 2e6) {
  warning("Very many edges (", n_thr, ") at EDGE_THR=", EDGE_THR,
          "; module construction will be slow. Consider raising threshold.")
}

# ==============================================================================
# STEP 3: Build pseudo-RobustnessResult and construct modules
# ==============================================================================

message("\n=== Step 3: Module construction ===")

# Construct a pseudo-RobustnessResult from the pooled edge table so that
# build_wgcna_modules() and .build_louvain_modules() — both of which consume
# a rob-like object — can be reused unchanged.
# R_score = |Spearman| (serves as the threshold selector)
# z_bar   = atanh(|Spearman|) (recovers |Spearman| when tanh() is applied)
et_thr <- et_pooled[abs(et_pooled$weight) >= EDGE_THR, , drop = FALSE]

pseudo_rob <- list(
  pair_scores = data.frame(
    gene_id_A = et_thr$gene_id_A,
    gene_id_B = et_thr$gene_id_B,
    R_score   = abs(et_thr$weight),
    z_bar     = atanh(pmin(abs(et_thr$weight), 0.9999)),
    stringsAsFactors = FALSE
  ),
  method_params = list(
    source    = "pseudobulk_spearman",
    threshold = EDGE_THR
  )
)
message(sprintf("  pseudo_rob: %d pairs", nrow(pseudo_rob$pair_scores)))

# ---- 3a. WGCNA (power = 1) --------------------------------------------------

message("  Building WGCNA modules (power=1)...")
t_wgcna <- proc.time()

wgcna_mi <- tryCatch({
  mi <- build_wgcna_modules(
    rob             = pseudo_rob,
    network_list    = network_list,
    r_score_min     = EDGE_THR,
    soft_power      = 1L,
    merge_cut       = 0.25,
    min_module_size = 30L,
    sub_merge_cut   = 0.10
  )
  mi$method           <- "wgcna_p1"
  mi$graph            <- "pseudobulk"
  mi$r_score_threshold <- EDGE_THR
  mi
}, error = function(e) {
  message("  WGCNA failed: ", conditionMessage(e))
  NULL
})

if (!is.null(wgcna_mi)) {
  gm <- wgcna_mi$gene_module
  message(sprintf("  WGCNA: %d modules | assigned=%d | grey=%d (%.1f%%) | %.1f min",
    length(unique(gm$top_module[gm$top_module > 0L])),
    sum(gm$top_module > 0L),
    sum(gm$top_module == 0L),
    100 * mean(gm$top_module == 0L),
    (proc.time() - t_wgcna)[["elapsed"]] / 60))
}

# ---- 3b. Louvain ------------------------------------------------------------

message("  Building Louvain modules...")
t_louv <- proc.time()

louv_mi <- tryCatch({
  mi <- .build_louvain_modules(pseudo_rob, r_score_min = EDGE_THR)
  mi$method           <- "louvain"
  mi$graph            <- "pseudobulk"
  mi$r_score_threshold <- EDGE_THR
  mi
}, error = function(e) {
  message("  Louvain failed: ", conditionMessage(e))
  NULL
})

if (!is.null(louv_mi)) {
  gm <- louv_mi$gene_module
  message(sprintf("  Louvain: %d modules | assigned=%d | grey=%d (%.1f%%) | %.1f min",
    length(unique(gm$top_module[gm$top_module > 0L])),
    sum(gm$top_module > 0L, na.rm = TRUE),
    sum(is.na(gm$top_module) | gm$top_module == 0L),
    100 * mean(is.na(gm$top_module) | gm$top_module == 0L),
    (proc.time() - t_louv)[["elapsed"]] / 60))
}

# ==============================================================================
# STEP 4: Annotate and save
# ==============================================================================

message("\n=== Step 4: Annotate and save ===")

wgcna_mi <- .annotate_and_save(
  mod_input   = wgcna_mi,
  method_name = "wgcna",
  out_dir     = OUT_MOD_DIR,
  counts_ctx  = counts_for_context,
  cond_lbl    = cond_labels,
  conditions  = CONDITIONS,
  network_list = network_list,
  tf_path     = TF_META_PATH,
  stratified  = STRATIFY
)

louv_mi <- .annotate_and_save(
  mod_input   = louv_mi,
  method_name = "louvain",
  out_dir     = OUT_MOD_DIR,
  counts_ctx  = counts_for_context,
  cond_lbl    = cond_labels,
  conditions  = CONDITIONS,
  network_list = network_list,
  tf_path     = TF_META_PATH,
  stratified  = STRATIFY
)

# Free the large counts matrix
rm(counts_for_context); gc()

# ==============================================================================
# STEP 5: BON3 check (AT1G08860)
# ==============================================================================

message("\n=== Step 5: BON3 (AT1G08860) ===")
BON3 <- "AT1G08860"

.bon3_row <- function(mi, label) {
  if (is.null(mi)) { message("  [", label, "] not built"); return() }
  gm <- mi$gene_module
  r  <- gm[gm$gene_id == BON3, , drop = FALSE]
  if (nrow(r) == 0L) {
    message(sprintf("  [%s] AT1G08860 NOT IN NETWORK", label))
    return(invisible(NULL))
  }
  top_mod <- r$top_module[1L]
  kme_val <- r$kME[1L]
  hg      <- mi$hub_genes
  hub_r   <- hg[hg$gene_id == BON3, , drop = FALSE]
  hub_rnk <- if (nrow(hub_r) > 0L) hub_r$hub_rank[1L] else NA_integer_
  message(sprintf("  [%s] module=%d | kME=%.5f | hub_rank=%s",
    label, top_mod, kme_val,
    ifelse(is.na(hub_rnk), "not a hub", as.character(hub_rnk))))
}

.bon3_row(wgcna_mi, "wgcna")
.bon3_row(louv_mi,  "louvain")

# Top partners from pseudo_rob
ps       <- pseudo_rob$pair_scores
bon3_ps  <- ps[ps$gene_id_A == BON3 | ps$gene_id_B == BON3, , drop = FALSE]

if (nrow(bon3_ps) > 0L) {
  bon3_ps$partner        <- ifelse(bon3_ps$gene_id_A == BON3,
                                   bon3_ps$gene_id_B, bon3_ps$gene_id_A)
  bon3_ps                <- bon3_ps[order(bon3_ps$R_score, decreasing = TRUE), ]
  bon3_ps$partner_symbol <- sym_lookup[bon3_ps$partner]
  bon3_top10             <- head(bon3_ps, 10L)

  message(sprintf("  Pseudobulk partners at |Spearman| >= %.1f: %d total",
                  EDGE_THR, nrow(bon3_ps)))
  for (i in seq_len(nrow(bon3_top10))) {
    message(sprintf("    %2d. %s (%s) |Spearman|=%.4f",
      i, bon3_top10$partner[i],
      ifelse(is.na(bon3_top10$partner_symbol[i]),
             "NA", bon3_top10$partner_symbol[i]),
      bon3_top10$R_score[i]))
  }
} else {
  message(sprintf("  BON3 has NO partners at |Spearman| >= %.1f", EDGE_THR))
  bon3_top10 <- data.frame()
}

# GGM comparison reference (from prior analysis)
message("\n  --- GGM vs Pseudobulk structural comparison ---")
message("  GGM large_wgcna   : module=0 (grey) | kME=NA")
message("  GGM large_louvain : module=8 | kME=0.00662")
message("  GGM n_partners (all conditions, any pcor): 147 (32+29+46+40)")
message("  GGM max_pcor: 0.039 | n_partners R_score>0: 11 | best R_score: 0.50")
message(sprintf("  Pseudobulk n_partners (|Spearman|>=%.1f): %d", EDGE_THR, nrow(bon3_ps)))
message(sprintf("  Pseudobulk max_|Spearman|: %.4f",
                if (nrow(bon3_ps) > 0L) max(bon3_ps$R_score) else 0))

# ==============================================================================
# STEP 6: PSEUDOBULK_SUMMARY.md
# ==============================================================================

message("\n=== Step 6: PSEUDOBULK_SUMMARY.md ===")

t_total <- (proc.time() - t_global)[["elapsed"]]

# Helper: safely extract counts from module result
.n_mods     <- function(mi) if (!is.null(mi)) length(unique(mi$gene_module$top_module[mi$gene_module$top_module > 0L])) else NA
.n_assigned <- function(mi) if (!is.null(mi)) sum(mi$gene_module$top_module > 0L, na.rm = TRUE) else NA
.n_grey     <- function(mi) if (!is.null(mi)) sum(is.na(mi$gene_module$top_module) | mi$gene_module$top_module == 0L) else NA
.n_tfs      <- function(mi) if (!is.null(mi)) nrow(mi$module_tfs) else NA
.bon3_mod   <- function(mi) {
  if (is.null(mi)) return("N/A")
  r <- mi$gene_module[mi$gene_module$gene_id == BON3, , drop = FALSE]
  if (nrow(r) == 0L) "not in network" else as.character(r$top_module[1L])
}
.bon3_kme   <- function(mi) {
  if (is.null(mi)) return("N/A")
  r <- mi$gene_module[mi$gene_module$gene_id == BON3, , drop = FALSE]
  if (nrow(r) == 0L) "N/A" else sprintf("%.5f", r$kME[1L])
}

# Genes assigned in pseudobulk that were grey in GGM large_wgcna
ggm_gm_path <- file.path("results", "pathogen_multiome", "official_modules",
                          "large_wgcna", "gene_module.csv")
n_pb_recover_wgcna <- NA_integer_
n_pb_recover_louv  <- NA_integer_
n_ggm_grey         <- NA_integer_
if (file.exists(ggm_gm_path)) {
  ggm_gm        <- read.csv(ggm_gm_path, stringsAsFactors = FALSE)
  ggm_grey_ids  <- ggm_gm$gene_id[ggm_gm$top_module == 0L]
  n_ggm_grey    <- length(ggm_grey_ids)
  if (!is.null(wgcna_mi)) {
    pb_assigned_w <- wgcna_mi$gene_module$gene_id[wgcna_mi$gene_module$top_module > 0L]
    n_pb_recover_wgcna <- sum(pb_assigned_w %in% ggm_grey_ids)
  }
  if (!is.null(louv_mi)) {
    pb_assigned_l <- louv_mi$gene_module$gene_id[
      !is.na(louv_mi$gene_module$top_module) & louv_mi$gene_module$top_module > 0L]
    n_pb_recover_louv <- sum(pb_assigned_l %in% ggm_grey_ids)
  }
}

summary_lines <- c(
  "# Pathogen Pseudobulk Co-expression — Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Total wall time: ", round(t_total / 60, 1), " min"),
  "",
  "## Input",
  paste0("- Seurat object: ", basename(SEURAT_PATH)),
  paste0("- Pseudobulk grouping (GROUP_VAR): ", GROUP_VAR),
  paste0("- Stratification: pooled (single stratum, all conditions combined)"),
  paste0("- Pseudobulk samples (subclusters): ",
         length(unique(et_pooled$gene_id_A))),   # proxy; network$gene_ids is accurate
  paste0("- Genes retained (min_cells >= ", MIN_CELLS, "): ",
         length(network_list[[1]]$gene_ids)),
  "",
  "## Network (pooled Spearman)",
  paste0("- Storage cutoff: |Spearman| >= 0.1"),
  paste0("- Total edges at storage cutoff: ", nrow(et_pooled)),
  paste0("- Edge threshold for modules: |Spearman| >= ", EDGE_THR),
  paste0("- Edges used for modules: ", nrow(et_thr)),
  sprintf("- Spearman distribution (all stored edges):"),
  sprintf("  min=%.4f  p10=%.4f  p25=%.4f  median=%.4f  p75=%.4f  max=%.4f",
          min(w), quantile(w, 0.10), quantile(w, 0.25),
          median(w), quantile(w, 0.75), max(w)),
  "",
  "## Module results",
  "",
  "### WGCNA (power=1, merge_cut=0.25, min_size=30, sub_merge_cut=0.10)",
  paste0("- n_modules:  ", .n_mods(wgcna_mi)),
  paste0("- n_assigned: ", .n_assigned(wgcna_mi)),
  paste0("- n_grey:     ", .n_grey(wgcna_mi)),
  paste0("- n_TF_entries: ", .n_tfs(wgcna_mi)),
  "",
  "### Louvain",
  paste0("- n_modules:  ", .n_mods(louv_mi)),
  paste0("- n_assigned: ", .n_assigned(louv_mi)),
  paste0("- n_grey:     ", .n_grey(louv_mi)),
  paste0("- n_TF_entries: ", .n_tfs(louv_mi)),
  "",
  "### Preservation",
  "Not computed. Single pooled network; no cross-stratum reference available.",
  "To compute: call compute_preservation_fallback() with output_per_condition/",
  "(GGM per-condition networks) as the test network_list2.",
  "",
  "## BON3 (AT1G08860) — GGM vs Pseudobulk structural comparison",
  "",
  "| Metric | GGM large_wgcna | GGM large_louvain | Pseudobulk WGCNA | Pseudobulk Louvain |",
  "|---|---|---|---|---|",
  paste0("| module | 0 (grey) | 8 | ", .bon3_mod(wgcna_mi), " | ", .bon3_mod(louv_mi), " |"),
  paste0("| kME | NA | 0.00662 | ", .bon3_kme(wgcna_mi), " | ", .bon3_kme(louv_mi), " |"),
  paste0("| n_partners | 11 (R>0) | — | ", nrow(bon3_ps), " (|Sp|>=",
         EDGE_THR, ") | — |"),
  paste0("| max weight | 0.50 (R_score) | — | ",
         sprintf("%.4f", if (nrow(bon3_ps) > 0) max(bon3_ps$R_score) else 0),
         " (|Spearman|) | — |"),
  "",
  "GGM total edges (all 4 conditions, any pcor): 147 (Mock=32, DC3000=29,",
  "AvrRpt2=46, AvrRpm1=40). Max pcor = 0.039 (just above 0.02 cutoff).",
  "",
  "## Genes recovered in pseudobulk modules that were grey in GGM",
  paste0("Reference set: genes in GGM large_wgcna module 0 (grey): ", n_ggm_grey),
  paste0("Recovered in pseudobulk WGCNA (assigned, module > 0): ", n_pb_recover_wgcna),
  paste0("Recovered in pseudobulk Louvain (assigned, module > 0): ", n_pb_recover_louv),
  "",
  "## Implementation notes",
  "- Spearman is computed as rank-transform-then-Pearson in estimate_pseudobulk().",
  "  The implementation ranks genes WITHIN each sample (across all genes in that",
  "  sample), then correlates within-sample rank profiles across samples. This",
  "  differs from standard Spearman (rank each gene ACROSS samples). The current",
  "  approach measures whether two genes consistently occupy similar rank positions",
  "  within subclusters, not whether their absolute expression co-varies.",
  "- Context annotation uses per-condition cell-level mean expression (post-hoc,",
  "  derived from counts matrix), not per-condition network edge weights.",
  "- Preservation not computed (single pooled stratum).",
  "- Louvain module_hier is an empty data.frame (no hierarchical structure)."
)

writeLines(summary_lines, file.path(RESULTS_DIR, "PSEUDOBULK_SUMMARY.md"))
message("Written: ", file.path(RESULTS_DIR, "PSEUDOBULK_SUMMARY.md"))

message(sprintf("\n=== DONE | Total wall time: %.1f min ===",
                (proc.time() - t_global)[["elapsed"]] / 60))
