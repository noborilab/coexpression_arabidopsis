## Official module construction — pseudobulk co-expression, pathogen_multiome
## Stage 3 threshold confirmed: global |r| = tanh(z_bar) >= 0.42
## Design: zscore_gene + Spearman + obs_subcluster (298 obs-points)
##
## Steps:
##   1.  Build edges_absr042.csv from pair_scores_full.csv
##   2.  Build sparse adjacency (Matrix::sparseMatrix)
##   3.  WGCNA pickSoftThreshold on obs-point expression matrix
##   4.  WGCNA blockwiseModules (90-min timeout; fallback power=6)
##   5.  Louvain modules (seed=98) + kME from expression
##   6.  GO enrichment (BP, BH p.adj < 0.05) for both
##   7.  Condition-pattern profiles for both
##   8.  Feature plots as PDFs
##   9.  BON3 / WRKY post-hoc sanity
##   10. Final report

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(igraph)
  library(jsonlite)
})

t_global <- proc.time()

# ==============================================================================
# PARAMETERS
# ==============================================================================

MIN_ABS_R        <- 0.42
SEED             <- 98L
MIN_MODULE_SIZE  <- 30L
MERGE_CUT        <- 0.25
WGCNA_POWER_VEC  <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20)
MAX_BLOCK_SIZE   <- 6000L
CONDITIONS       <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")
BON3_ID          <- "AT1G08860"
DATASET_ID       <- "pathogen_multiome"

RESULTS_DIR  <- file.path("results", DATASET_ID)
PB_DIR       <- file.path(RESULTS_DIR, "pseudobulk_zscore_spearman")
OUT_DIR      <- file.path(PB_DIR, "modules_official")
STAGE3_DIR   <- file.path(RESULTS_DIR, "stage3_threshold_sweep")
WGCNA_DIR    <- file.path(OUT_DIR, "wgcna")
LOUVAIN_DIR  <- file.path(OUT_DIR, "louvain")
PLOTS_DIR    <- file.path(OUT_DIR, "plots")

PAIR_SCORES_PATH <- file.path(PB_DIR, "pair_scores_full.csv")
OBS_CACHE_PATH   <- file.path(STAGE3_DIR, "obs_normalized_cache.rds")
SYMBOL_MAP_PATH  <- file.path(RESULTS_DIR, "symbol_map.csv")
WRKY_PATH        <- file.path(RESULTS_DIR, "geneset_lookups", "WRKY_GGM_vs_PB.csv")
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

for (d in c(OUT_DIR, WGCNA_DIR, LOUVAIN_DIR, PLOTS_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

symbol_map <- read.csv(SYMBOL_MAP_PATH, stringsAsFactors = FALSE)
sym_lookup  <- setNames(symbol_map$gene_symbol, symbol_map$gene_id)
message("symbol_map: ", nrow(symbol_map), " entries")

# ==============================================================================
# STEP 1: Build edges_absr042.csv
# ==============================================================================

message("\n==== STEP 1: Build |r|>=", MIN_ABS_R, " edge list ====")
t1 <- proc.time()

EDGES_PATH <- file.path(OUT_DIR, "edges_absr042.csv")

ps_cols <- c("gene_id_A", "gene_id_B", "z_bar", "R_score",
             "I_Mock", "I_DC3000", "I_AvrRpt2", "I_AvrRpm1")
message("  Reading pair_scores_full.csv (", round(file.info(PAIR_SCORES_PATH)$size / 1e9, 1), " GB)...")
ps_all <- data.table::fread(PAIR_SCORES_PATH, nThread = 1L,
                            select = ps_cols, data.table = FALSE)
message("  Loaded: ", nrow(ps_all), " pairs")

ps_all$mean_abs_r <- abs(tanh(ps_all$z_bar))

keep_idx <- which(!is.na(ps_all$mean_abs_r) & ps_all$mean_abs_r >= MIN_ABS_R)
ps042 <- ps_all[keep_idx, , drop = FALSE]
rm(ps_all, keep_idx); invisible(gc())

# Enforce alphabetical order gene_id_A < gene_id_B
swap_idx <- which(ps042$gene_id_A > ps042$gene_id_B)
if (length(swap_idx) > 0L) {
  tmp                      <- ps042$gene_id_A[swap_idx]
  ps042$gene_id_A[swap_idx] <- ps042$gene_id_B[swap_idx]
  ps042$gene_id_B[swap_idx] <- tmp
  rm(tmp)
}
pair_key <- paste(ps042$gene_id_A, ps042$gene_id_B, sep = "__")
ps042 <- ps042[!duplicated(pair_key), , drop = FALSE]
rm(pair_key)

n_pairs <- nrow(ps042)
visible_genes <- sort(unique(c(ps042$gene_id_A, ps042$gene_id_B)))
n_genes       <- length(visible_genes)

edges_df <- data.frame(
  gene_id_A  = ps042$gene_id_A,
  gene_id_B  = ps042$gene_id_B,
  mean_abs_r = round(ps042$mean_abs_r, 8),
  stringsAsFactors = FALSE
)
write.csv(edges_df, EDGES_PATH, row.names = FALSE)

cat(sprintf("STEP 1 complete — n_pairs=%d (expected~751,959)  n_genes=%d (expected~5,450)\n",
            n_pairs, n_genes))

# ==============================================================================
# STEP 2: Build sparse adjacency matrix
# ==============================================================================

message("\n==== STEP 2: Build sparse adjacency matrix ====")
t2 <- proc.time()

gene_idx <- setNames(seq_along(visible_genes), visible_genes)
i_idx    <- gene_idx[ps042$gene_id_A]
j_idx    <- gene_idx[ps042$gene_id_B]

A_sparse <- Matrix::sparseMatrix(
  i    = c(i_idx, j_idx),
  j    = c(j_idx, i_idx),
  x    = c(edges_df$mean_abs_r, edges_df$mean_abs_r),
  dims = c(n_genes, n_genes),
  dimnames = list(visible_genes, visible_genes)
)
rm(i_idx, j_idx)

cat(sprintf("STEP 2 complete — sparse adjacency %dx%d  nnz=%d\n",
            nrow(A_sparse), ncol(A_sparse), Matrix::nnzero(A_sparse)))

# ==============================================================================
# STEP 3: WGCNA pickSoftThreshold
# ==============================================================================

message("\n==== STEP 3: WGCNA pickSoftThreshold ====")
t3 <- proc.time()

if (!requireNamespace("WGCNA", quietly = TRUE))
  stop("WGCNA not installed. Run: BiocManager::install('WGCNA')")

suppressPackageStartupMessages(library(WGCNA))
WGCNA::allowWGCNAThreads()

message("  Loading obs_normalized_cache.rds...")
obs_cache  <- readRDS(OBS_CACHE_PATH)
obs_mat    <- obs_cache$matrix                         # 11010 × 298 (genes × obs-pts)
vis_in_cache <- intersect(visible_genes, rownames(obs_mat))
message("  Visible genes in obs cache: ", length(vis_in_cache), " / ", n_genes)

obs_masked <- obs_mat[vis_in_cache, , drop = FALSE]   # n_vis × 298
rm(obs_mat, obs_cache)
datExpr    <- t(obs_masked)                            # 298 × n_vis (samples × genes)
message(sprintf("  datExpr: %d obs-points × %d genes", nrow(datExpr), ncol(datExpr)))

message("  Running pickSoftThreshold...")
set.seed(SEED)
sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = WGCNA_POWER_VEC,
  networkType = "unsigned",
  verbose     = 5
)

fit_tbl <- sft$fitIndices
message("  Fit table:")
print(fit_tbl)

r2_col   <- "SFT.R.sq"
good_pow <- fit_tbl[!is.na(fit_tbl[[r2_col]]) & fit_tbl[[r2_col]] >= 0.80, ]
if (nrow(good_pow) > 0L) {
  selected_power <- good_pow$Power[1L]
  r2_selected    <- good_pow[[r2_col]][1L]
  r2_threshold_met <- TRUE
  message(sprintf("  Selected power=%d  R²=%.3f (first >= 0.80)", selected_power, r2_selected))
} else {
  r2_vals  <- fit_tbl[[r2_col]]
  r2_vals[is.na(r2_vals)] <- 0
  elbow_i  <- which.max(diff(r2_vals))
  selected_power    <- fit_tbl$Power[elbow_i]
  r2_selected       <- r2_vals[elbow_i]
  r2_threshold_met  <- FALSE
  warning_msg <- sprintf("WARNING: scale-free R2 did not reach 0.80; using elbow power=%d (R²=%.3f)",
                         selected_power, r2_selected)
  message("  ", warning_msg)
  cat(warning_msg, "\n")
}

# Save plot
plot_path <- file.path(OUT_DIR, "pickSoftThreshold_plot.pdf")
pdf(plot_path, width = 10, height = 5)
par(mfrow = c(1, 2))
plot(fit_tbl$Power, fit_tbl[[r2_col]], type = "b", pch = 19,
     xlab = "Soft Power", ylab = "Scale Free Topology R²",
     main = "Scale-free R² vs Power")
abline(h = 0.80, col = "red", lty = 2)
abline(v = selected_power, col = "blue", lty = 2)
text(selected_power, max(r2_selected + 0.03, 0.1),
     paste0("p=", selected_power), col = "blue", cex = 0.8)

plot(fit_tbl$Power, fit_tbl$mean.k., type = "b", pch = 19,
     xlab = "Soft Power", ylab = "Mean Connectivity",
     main = "Mean Connectivity vs Power")
abline(v = selected_power, col = "blue", lty = 2)
dev.off()
message("  Plot saved: ", plot_path)

write.csv(fit_tbl, file.path(OUT_DIR, "pickSoftThreshold_fitIndices.csv"), row.names = FALSE)

cat(sprintf("STEP 3 complete — power=%d  R²=%.3f  R2>=0.80=%s\n",
            selected_power, r2_selected, if (r2_threshold_met) "YES" else "NO (elbow)"))

# ==============================================================================
# STEP 4: WGCNA blockwiseModules
# ==============================================================================

message("\n==== STEP 4: WGCNA blockwiseModules (power=", selected_power, ") ====")
t4 <- proc.time()

wgcna_net       <- NULL
wgcna_ok        <- FALSE
wgcna_power     <- selected_power
wgcna_labels_out <- NULL
wgcna_n_mods    <- NA_integer_
wgcna_grey_rate <- NA_real_
sizes_w         <- integer(0)

run_blockwise <- function(pow) {
  set.seed(SEED)
  WGCNA::blockwiseModules(
    datExpr,
    power             = pow,
    networkType       = "unsigned",
    TOMType           = "unsigned",
    minModuleSize     = MIN_MODULE_SIZE,
    mergeCutHeight    = MERGE_CUT,
    numericLabels     = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs          = FALSE,
    verbose           = 3,
    maxBlockSize      = MAX_BLOCK_SIZE
  )
}

# 90-minute time box
message("  Running blockwiseModules (90-min limit)...")
setTimeLimit(elapsed = 5400)
wgcna_net <- tryCatch({
  r <- run_blockwise(selected_power)
  setTimeLimit()
  wgcna_ok <- TRUE
  r
}, error = function(e) {
  setTimeLimit()
  msg <- conditionMessage(e)
  if (grepl("elapsed time limit", msg, ignore.case = TRUE)) {
    message("  blockwiseModules timeout — retrying with fallback power=6")
    wgcna_power <<- 6L
    setTimeLimit(elapsed = 5400)
    r2 <- tryCatch({
      res <- run_blockwise(6L)
      setTimeLimit()
      wgcna_ok <<- TRUE
      res
    }, error = function(e2) {
      setTimeLimit()
      message("  Fallback also failed: ", conditionMessage(e2))
      NULL
    })
    r2
  } else {
    message("  blockwiseModules error: ", msg)
    NULL
  }
})
setTimeLimit()

if (!is.null(wgcna_net) && wgcna_ok) {
  net_colors  <- wgcna_net$colors     # named integer vector, 0 = grey
  wgcna_genes <- names(net_colors)

  # Module eigengenes and signed kME
  wgcna_MEs <- WGCNA::moduleEigengenes(
    expr    = datExpr[, wgcna_genes, drop = FALSE],
    colors  = net_colors,
    verbose = 0
  )$eigengenes

  kme_mat <- WGCNA::signedKME(
    datExpr[, wgcna_genes, drop = FALSE],
    wgcna_MEs,
    outputColumnName = "kME"
  )

  # Extract per-gene kME for assigned module
  # signedKME with outputColumnName="kME" produces "kME0","kME1",... (NOT "kMEME0")
  kme_vec <- vapply(seq_along(wgcna_genes), function(i) {
    mod <- net_colors[i]
    if (mod == 0L) return(NA_real_)
    col_nm <- paste0("kME", mod)
    if (col_nm %in% colnames(kme_mat)) kme_mat[i, col_nm] else NA_real_
  }, numeric(1L))

  wgcna_labels_out <- data.frame(
    gene_id = wgcna_genes,
    module  = as.integer(net_colors),
    kME     = kme_vec,
    stringsAsFactors = FALSE
  )

  mod_counts_w <- table(wgcna_labels_out$module[wgcna_labels_out$module > 0L])
  wgcna_n_mods <- length(mod_counts_w)
  n_grey_w     <- sum(wgcna_labels_out$module == 0L)
  wgcna_grey_rate <- n_grey_w / nrow(wgcna_labels_out)
  sizes_w      <- sort(as.integer(mod_counts_w))

  get_top5 <- function(m) {
    sub <- wgcna_labels_out[wgcna_labels_out$module == m & !is.na(wgcna_labels_out$kME), ]
    sub <- sub[order(sub$kME, decreasing = TRUE), ]
    paste(head(sub$gene_id, 5L), collapse = ";")
  }
  wgcna_summary <- data.frame(
    module         = as.integer(names(mod_counts_w)),
    n_genes        = as.integer(mod_counts_w),
    top5_kME_genes = vapply(as.integer(names(mod_counts_w)), get_top5, character(1L)),
    stringsAsFactors = FALSE
  )

  write.csv(wgcna_labels_out, file.path(WGCNA_DIR, "module_membership.csv"), row.names = FALSE)
  write.csv(wgcna_summary,    file.path(WGCNA_DIR, "module_summary.csv"),    row.names = FALSE)
  saveRDS(wgcna_net, file.path(WGCNA_DIR, "blockwiseModules_net.rds"))

  wgcna_params <- list(
    power          = wgcna_power,
    networkType    = "unsigned",
    minModuleSize  = MIN_MODULE_SIZE,
    mergeCutHeight = MERGE_CUT,
    n_modules      = wgcna_n_mods,
    grey_rate      = round(wgcna_grey_rate, 4),
    timestamp      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  writeLines(jsonlite::toJSON(wgcna_params, pretty = TRUE, auto_unbox = TRUE),
             file.path(WGCNA_DIR, "wgcna_params.json"))

  message(sprintf(
    "  WGCNA: n_modules=%d  grey=%d/%d (%.1f%%)  size min/med/max=%d/%.0f/%d  power=%d",
    wgcna_n_mods, n_grey_w, nrow(wgcna_labels_out), 100 * wgcna_grey_rate,
    min(sizes_w), median(sizes_w), max(sizes_w), wgcna_power
  ))
}

cat(sprintf("STEP 4 complete — WGCNA: n_modules=%s  grey_rate=%s  power=%d\n",
            if (is.na(wgcna_n_mods)) "FAILED" else as.character(wgcna_n_mods),
            if (is.na(wgcna_grey_rate)) "NA" else sprintf("%.1f%%", 100 * wgcna_grey_rate),
            wgcna_power))

# ==============================================================================
# STEP 5: Louvain module construction
# ==============================================================================

message("\n==== STEP 5: Louvain modules (seed=", SEED, ") ====")
t5 <- proc.time()

g_louvain <- igraph::graph_from_data_frame(
  d = data.frame(from = edges_df$gene_id_A, to = edges_df$gene_id_B,
                 weight = edges_df$mean_abs_r, stringsAsFactors = FALSE),
  directed = FALSE,
  vertices = data.frame(name = visible_genes, stringsAsFactors = FALSE)
)

set.seed(SEED)
cl <- igraph::cluster_louvain(g_louvain, weights = igraph::E(g_louvain)$weight)
rm(g_louvain)

memb   <- igraph::membership(cl)
lv_genes <- names(memb)
lv_lbl   <- as.integer(memb)

# Small communities → module 0 (grey)
comm_sz <- table(lv_lbl)
small_c <- as.integer(names(comm_sz[comm_sz < MIN_MODULE_SIZE]))
lv_lbl[lv_lbl %in% small_c] <- 0L

# Relabel surviving modules as consecutive integers
live_c  <- sort(unique(lv_lbl[lv_lbl > 0L]))
relabel <- setNames(seq_along(live_c), as.character(live_c))
lv_lbl  <- ifelse(lv_lbl == 0L, 0L, as.integer(relabel[as.character(lv_lbl)]))
names(lv_lbl) <- lv_genes

# kME: correlation of gene's expression profile with its module centroid
louvain_kme  <- rep(NA_real_, length(lv_genes))
names(louvain_kme) <- lv_genes
vis_in_cache_set <- vis_in_cache   # already computed in step 3

for (m in sort(unique(lv_lbl[lv_lbl > 0L]))) {
  m_genes_vis <- intersect(lv_genes[!is.na(lv_lbl) & lv_lbl == m], vis_in_cache_set)
  if (length(m_genes_vis) < 2L) next
  mod_mean <- colMeans(obs_masked[m_genes_vis, , drop = FALSE])
  for (g in m_genes_vis)
    louvain_kme[g] <- cor(obs_masked[g, ], mod_mean)
}

louvain_labels_out <- data.frame(
  gene_id = lv_genes,
  module  = as.integer(lv_lbl),
  kME     = louvain_kme[lv_genes],
  stringsAsFactors = FALSE
)

mod_counts_l <- table(louvain_labels_out$module[louvain_labels_out$module > 0L])
lv_n_mods    <- length(mod_counts_l)
n_grey_lv    <- sum(louvain_labels_out$module == 0L)
lv_grey_rate <- n_grey_lv / nrow(louvain_labels_out)
sizes_lv     <- sort(as.integer(mod_counts_l))

get_top5_lv <- function(m) {
  sub <- louvain_labels_out[louvain_labels_out$module == m & !is.na(louvain_labels_out$kME), ]
  sub <- sub[order(sub$kME, decreasing = TRUE), ]
  paste(head(sub$gene_id, 5L), collapse = ";")
}
louvain_summary <- data.frame(
  module         = as.integer(names(mod_counts_l)),
  n_genes        = as.integer(mod_counts_l),
  top5_kME_genes = vapply(as.integer(names(mod_counts_l)), get_top5_lv, character(1L)),
  stringsAsFactors = FALSE
)

write.csv(louvain_labels_out, file.path(LOUVAIN_DIR, "module_membership.csv"), row.names = FALSE)
write.csv(louvain_summary,    file.path(LOUVAIN_DIR, "module_summary.csv"),    row.names = FALSE)

louvain_params <- list(
  seed       = SEED,
  weight_col = "mean_abs_r",
  n_modules  = lv_n_mods,
  grey_rate  = round(lv_grey_rate, 4),
  timestamp  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)
writeLines(jsonlite::toJSON(louvain_params, pretty = TRUE, auto_unbox = TRUE),
           file.path(LOUVAIN_DIR, "louvain_params.json"))

message(sprintf("  Louvain: n_modules=%d  grey=%d/%d (%.1f%%)  size min/med/max=%d/%.0f/%d",
                lv_n_mods, n_grey_lv, nrow(louvain_labels_out), 100 * lv_grey_rate,
                min(sizes_lv), median(sizes_lv), max(sizes_lv)))

cat(sprintf("STEP 5 complete — Louvain: n_modules=%d  grey_rate=%.1f%%  size min/med/max=%d/%.0f/%d\n",
            lv_n_mods, 100 * lv_grey_rate,
            min(sizes_lv), median(sizes_lv), max(sizes_lv)))

# ==============================================================================
# STEP 6: GO enrichment
# ==============================================================================

message("\n==== STEP 6: GO enrichment ====")
t6 <- proc.time()

n_go_wgcna   <- NA_integer_
n_go_louvain <- NA_integer_
go_ok        <- requireNamespace("clusterProfiler", quietly = TRUE) &&
                requireNamespace("org.At.tair.db",  quietly = TRUE)

if (!go_ok) {
  message("  clusterProfiler or org.At.tair.db not available — skipping GO")
} else {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.At.tair.db)
  })

  run_go <- function(membership_df, out_dir) {
    bg_genes <- membership_df$gene_id
    mods     <- sort(unique(membership_df$module))
    mods     <- mods[!is.na(mods) & mods > 0L]
    go_rows  <- lapply(mods, function(m) {
      genes_m <- membership_df$gene_id[membership_df$module == m]
      tryCatch({
        res <- clusterProfiler::enrichGO(
          gene          = genes_m,
          universe      = bg_genes,
          OrgDb         = org.At.tair.db,
          keyType       = "TAIR",
          ont           = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.05,
          minGSSize     = 10,
          readable      = FALSE
        )
        if (is.null(res) || nrow(res@result) == 0L) return(NULL)
        df <- as.data.frame(res)
        data.frame(
          module      = m,
          GO_ID       = df$ID,
          Description = df$Description,
          p.adjust    = df$p.adjust,
          gene_ratio  = df$GeneRatio,
          bg_ratio    = df$BgRatio,
          gene_ids    = df$geneID,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        message("    GO failed module ", m, ": ", conditionMessage(e)); NULL
      })
    })
    go_df <- do.call(rbind, Filter(Negate(is.null), go_rows))
    if (is.null(go_df) || nrow(go_df) == 0L)
      go_df <- data.frame(module = integer(), GO_ID = character(),
                          Description = character(), p.adjust = numeric(),
                          gene_ratio = character(), bg_ratio = character(),
                          gene_ids = character(), stringsAsFactors = FALSE)
    write.csv(go_df, file.path(out_dir, "go_enrichment.csv"), row.names = FALSE)
    message("    Saved ", nrow(go_df), " terms → ", file.path(out_dir, "go_enrichment.csv"))
    go_df
  }

  if (!is.null(wgcna_labels_out)) {
    message("  GO for WGCNA...")
    go_wgcna     <- run_go(wgcna_labels_out, WGCNA_DIR)
    n_go_wgcna   <- length(unique(go_wgcna$module[go_wgcna$p.adjust < 0.05]))
  }
  message("  GO for Louvain...")
  go_louvain    <- run_go(louvain_labels_out, LOUVAIN_DIR)
  n_go_louvain  <- length(unique(go_louvain$module[go_louvain$p.adjust < 0.05]))
}

cat(sprintf("STEP 6 complete — GO: WGCNA %s mods  Louvain %s mods with ≥1 term\n",
            if (is.na(n_go_wgcna)) "SKIP" else as.character(n_go_wgcna),
            if (is.na(n_go_louvain)) "SKIP" else as.character(n_go_louvain)))

# ==============================================================================
# STEP 7: Condition-pattern profiles
# ==============================================================================

message("\n==== STEP 7: Condition-pattern profiles ====")
t7 <- proc.time()

# Compute pattern from I_ columns (avoids reading 9.8 GB pair_condition_patterns.csv)
I_cols      <- paste0("I_", CONDITIONS)
ps042$pattern <- apply(ps042[, I_cols, drop = FALSE], 1, function(x)
  paste(as.integer(x), collapse = ""))
ps042$n_cond  <- rowSums(ps042[, I_cols, drop = FALSE], na.rm = TRUE)

build_cond_profiles <- function(membership_df, edges_ps, out_dir) {
  gene_mod <- setNames(membership_df$module, membership_df$gene_id)
  mod_A    <- gene_mod[edges_ps$gene_id_A]
  mod_B    <- gene_mod[edges_ps$gene_id_B]
  same_mod <- !is.na(mod_A) & !is.na(mod_B) & mod_A == mod_B & mod_A > 0L
  ps_mod   <- edges_ps[same_mod, , drop = FALSE]
  ps_mod$module <- mod_A[same_mod]

  mods <- sort(unique(ps_mod$module))
  rows <- lapply(mods, function(m) {
    sub     <- ps_mod[ps_mod$module == m, , drop = FALSE]
    pat_tab <- table(sub$pattern)
    dom_pat <- names(which.max(pat_tab))
    pat_frac <- round(as.numeric(pat_tab) / sum(pat_tab), 4)
    names(pat_frac) <- paste0("frac_", names(pat_tab))
    # Use list() to preserve types (c() with character coerces numerics)
    row_base <- list(
      module           = m,
      n_pairs          = nrow(sub),
      dominant_pattern = dom_pat,
      mean_r_score     = round(mean(sub$R_score, na.rm = TRUE), 4),
      mean_n_cond      = round(mean(sub$n_cond,  na.rm = TRUE), 3)
    )
    c(row_base, as.list(pat_frac))
  })

  # Build output data frame robustly (rows may have different pattern columns)
  all_nms <- unique(unlist(lapply(rows, names)))
  out_list <- lapply(rows, function(r) {
    out <- as.list(r)
    miss <- setdiff(all_nms, names(out))
    for (nm in miss) out[[nm]] <- NA
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  out_df <- do.call(rbind, out_list)

  # Coerce numeric columns
  num_nms <- setdiff(names(out_df), "dominant_pattern")
  for (nm in num_nms)
    out_df[[nm]] <- suppressWarnings(as.numeric(out_df[[nm]]))

  write.csv(out_df, file.path(out_dir, "module_condition_patterns.csv"), row.names = FALSE)
  message("    Saved ", nrow(out_df), " rows → ",
          file.path(out_dir, "module_condition_patterns.csv"))
  invisible(out_df)
}

cp_wgcna   <- if (!is.null(wgcna_labels_out))
  build_cond_profiles(wgcna_labels_out,  ps042, WGCNA_DIR) else NULL
cp_louvain <- build_cond_profiles(louvain_labels_out, ps042, LOUVAIN_DIR)

cat(sprintf("STEP 7 complete — condition profiles: WGCNA %s rows  Louvain %s rows\n",
            if (is.null(cp_wgcna)) "SKIP" else as.character(nrow(cp_wgcna)),
            as.character(nrow(cp_louvain))))

# ==============================================================================
# STEP 8: Feature plots
# ==============================================================================

message("\n==== STEP 8: Feature plots ====")
t8 <- proc.time()

PURPLE_COLS <- c("lightgray", "#BFD3E6", "#9EBCDA", "#8C96C6",
                 "#8C6BB1", "#88419D", "#810F7C", "#4D004B")
TOP_N <- 4L

make_feature_plots <- function(membership_df, summary_df, set_label, pdf_path) {
  if (nrow(summary_df) == 0L) {
    message("  No modules for ", set_label, " — skipping feature plots")
    return(invisible(NULL))
  }
  if (!file.exists(SEURAT_PATH)) {
    message("  Seurat object not found — skipping feature plots for ", set_label)
    return(invisible(NULL))
  }
  suppressPackageStartupMessages({
    library(Seurat); library(ggplot2); library(patchwork)
  })
  message("  Loading Seurat object for ", set_label, "...")
  sobj <- tryCatch(readRDS(SEURAT_PATH), error = function(e) {
    message("  Seurat load failed: ", conditionMessage(e)); NULL
  })
  if (is.null(sobj)) return(invisible(NULL))

  s2   <- gsub("_rep[12]$", "", sobj$sample)
  s2   <- gsub("_(04|06|09|24)h$", "", s2)
  s2   <- gsub("^00_Mock$", "Mock", s2)
  sobj$condition <- factor(s2, levels = CONDITIONS)
  DefaultAssay(sobj) <- "RNA"
  RNA_GENES <- rownames(sobj[["RNA"]])

  .resolve <- function(gid) {
    sym   <- sym_lookup[gid]
    cands <- c(if (!is.na(sym)) sym else character(0), gid)
    for (cand in cands) if (cand %in% RNA_GENES) return(cand)
    NA_character_
  }
  .lbl <- function(gid) {
    sym <- sym_lookup[gid]
    if (!is.na(sym) && nchar(trimws(sym)) > 0L) paste0(gid, " / ", sym) else gid
  }

  all_mods <- sort(unique(summary_df$module))
  message(sprintf("  %s: %d modules → %s", set_label, length(all_mods), pdf_path))
  pdf(pdf_path, width = length(CONDITIONS) * 5,
      height = (TOP_N + 1L) * 5, onefile = TRUE)

  for (mid in all_mods) {
    sub      <- membership_df[membership_df$module == mid & !is.na(membership_df$kME), ]
    sub      <- head(sub[order(sub$kME, decreasing = TRUE), ], TOP_N)
    hub_gids <- sub$gene_id
    feats    <- vapply(hub_gids, .resolve, character(1L))
    ok       <- !is.na(feats)
    if (!any(ok)) next

    all_gids <- membership_df$gene_id[membership_df$module == mid]
    all_rna  <- na.omit(vapply(all_gids, .resolve, character(1L)))
    if (length(all_rna) == 0L) next

    tryCatch({
      tmp        <- AddModuleScore(sobj, features = list(all_rna), name = "TmpMS_")
      mod_scores <- tmp@meta.data[["TmpMS_1"]]
      rm(tmp)

      umap_df <- as.data.frame(Embeddings(sobj, reduction = "umap"))
      colnames(umap_df) <- c("UMAP_1", "UMAP_2")
      umap_df$condition    <- sobj$condition
      umap_df$module_score <- mod_scores

      gene_rows <- lapply(which(ok), function(i) {
        fp <- FeaturePlot(sobj, features = feats[i], split.by = "condition",
                          reduction = "umap", pt.size = 0.5, order = TRUE,
                          max.cutoff = "q99", min.cutoff = "q1")
        fp <- fp & scale_colour_gradientn(colours = PURPLE_COLS, na.value = "lightgray") & NoAxes()
        for (k in seq_along(fp)) {
          ttl <- if (k == 1L) paste0(.lbl(hub_gids[i]), "\n(", CONDITIONS[k], ")") else CONDITIONS[k]
          fp[[k]] <- fp[[k]] + labs(title = ttl) +
            theme(plot.title = element_text(size = 9, face = "bold"))
        }
        fp
      })

      score_panels <- lapply(seq_along(CONDITIONS), function(k) {
        df_k <- umap_df[umap_df$condition == CONDITIONS[k], ]
        df_k <- df_k[order(df_k$module_score, na.last = FALSE), ]
        ggplot(df_k, aes(x = UMAP_1, y = UMAP_2, color = module_score)) +
          geom_point(size = 0.5) +
          scale_color_gradientn(colours = PURPLE_COLS, na.value = "lightgray", name = "score") +
          labs(title = if (k == 1L) paste0("M", mid, " score\n(", CONDITIONS[k], ")") else CONDITIONS[k]) +
          NoAxes() +
          theme(plot.background = element_blank(), panel.background = element_blank(),
                plot.title = element_text(size = 9, face = "bold"))
      })

      print(wrap_plots(c(gene_rows, list(wrap_plots(score_panels, nrow = 1L))), ncol = 1L))
    }, error = function(e)
      message(sprintf("  ERROR M%s: %s", mid, conditionMessage(e)))
    )
  }
  dev.off()
  rm(sobj); invisible(gc())
  message("  Saved: ", pdf_path)
}

make_feature_plots(
  membership_df = if (!is.null(wgcna_labels_out)) wgcna_labels_out
                  else data.frame(gene_id=character(), module=integer(), kME=numeric()),
  summary_df    = if (!is.null(wgcna_labels_out)) wgcna_summary
                  else data.frame(module=integer()),
  set_label     = "wgcna",
  pdf_path      = file.path(PLOTS_DIR, "featureplots_wgcna.pdf")
)

make_feature_plots(
  membership_df = louvain_labels_out,
  summary_df    = louvain_summary,
  set_label     = "louvain",
  pdf_path      = file.path(PLOTS_DIR, "featureplots_louvain.pdf")
)

cat(sprintf("STEP 8 complete — plots → %s\n", PLOTS_DIR))

# ==============================================================================
# STEP 9: BON3 / WRKY post-hoc sanity
# ==============================================================================

message("\n==== STEP 9: BON3 / WRKY post-hoc sanity ====")
t9 <- proc.time()

wrky_ids <- character(0)
if (file.exists(WRKY_PATH)) {
  wrky_df  <- read.csv(WRKY_PATH, stringsAsFactors = FALSE)
  wrky_ids <- wrky_df$gene_id
  message("  WRKY IDs: ", length(wrky_ids))
} else {
  message("  WRKY_PATH not found: ", WRKY_PATH)
}

bon3_042 <- ps042[ps042$gene_id_A == BON3_ID | ps042$gene_id_B == BON3_ID, , drop = FALSE]
bon3_042$partner <- ifelse(bon3_042$gene_id_A == BON3_ID, bon3_042$gene_id_B, bon3_042$gene_id_A)
bon3_042 <- bon3_042[order(bon3_042$mean_abs_r, decreasing = TRUE), ]

.report_gene <- function(gid, mem, set_lbl) {
  r <- mem[mem$gene_id == gid, , drop = FALSE]
  if (nrow(r) == 0L) {
    cat(sprintf("    [%s] not in network\n", set_lbl))
  } else {
    mod <- r$module[1L]
    kme <- r$kME[1L]
    cat(sprintf("    [%s] module=%s  kME=%s\n",
                set_lbl,
                if (is.na(mod) || mod == 0L) "grey(0)" else as.character(mod),
                if (is.na(kme)) "NA" else sprintf("%.5f", kme)))
  }
}

empty_mem <- data.frame(gene_id = character(), module = integer(), kME = numeric(),
                        stringsAsFactors = FALSE)

cat("\n===== BON3 / WRKY Post-hoc Sanity (NOT selection input) =====\n\n")
cat(sprintf("BON3 (%s):\n", BON3_ID))
cat(sprintf("  n_partners at |r|>=%.2f: %d\n", MIN_ABS_R, nrow(bon3_042)))
.report_gene(BON3_ID, if (!is.null(wgcna_labels_out)) wgcna_labels_out else empty_mem, "WGCNA")
.report_gene(BON3_ID, louvain_labels_out, "Louvain")

if (nrow(bon3_042) > 0L) {
  top10 <- head(bon3_042, 10L)
  cat(sprintf("\n  BON3 top-%d partners (by |r|):\n", nrow(top10)))
  for (i in seq_len(nrow(top10))) {
    sym_s <- sym_lookup[top10$partner[i]]
    cat(sprintf("    %2d. %s (%s)  |r|=%.4f  R_score=%.0f\n",
                i, top10$partner[i], ifelse(is.na(sym_s), "NA", sym_s),
                top10$mean_abs_r[i], top10$R_score[i]))
  }
}

.report_wrky <- function(mem, set_lbl) {
  if (length(wrky_ids) == 0L || nrow(mem) == 0L) return(invisible(NULL))
  gm_w  <- mem[mem$gene_id %in% wrky_ids, , drop = FALSE]
  n_net <- nrow(gm_w)
  n_asn <- sum(!is.na(gm_w$module) & gm_w$module > 0L)
  kme_a <- gm_w$kME[!is.na(gm_w$kME) & !is.na(gm_w$module) & gm_w$module > 0L]
  kq    <- if (length(kme_a) > 0L) quantile(kme_a, c(0, .25, .5, .75, 1)) else rep(NA_real_, 5)
  cat(sprintf(
    "\n  [%s] n_in_network=%d  n_assigned=%d  grey=%d\n  kME: min=%.3f Q1=%.3f med=%.3f Q3=%.3f max=%.3f\n",
    set_lbl, n_net, n_asn, n_net - n_asn,
    kq[1], kq[2], kq[3], kq[4], kq[5]
  ))
  top5w <- gm_w[!is.na(gm_w$kME) & !is.na(gm_w$module) & gm_w$module > 0L, ]
  top5w <- head(top5w[order(top5w$kME, decreasing = TRUE), ], 5L)
  if (nrow(top5w) > 0L) {
    cat(sprintf("  Top-5 WRKY by kME [%s]:\n", set_lbl))
    for (i in seq_len(nrow(top5w))) {
      sym_w <- sym_lookup[top5w$gene_id[i]]
      cat(sprintf("    %s (%s)  module=%d  kME=%.4f\n",
                  top5w$gene_id[i], ifelse(is.na(sym_w), "NA", sym_w),
                  top5w$module[i], top5w$kME[i]))
    }
  }
}

cat(sprintf("\nWRKY family (%d AT-IDs):\n", length(wrky_ids)))
.report_wrky(if (!is.null(wgcna_labels_out)) wgcna_labels_out else empty_mem, "WGCNA")
.report_wrky(louvain_labels_out, "Louvain")

cat(sprintf("\nSTEP 9 complete — BON3 partners=%d  WRKY_IDs=%d\n",
            nrow(bon3_042), length(wrky_ids)))

# ==============================================================================
# STEP 10: Final report
# ==============================================================================

t_total <- (proc.time() - t_global)[["elapsed"]]

cat("\n\n========== FINAL REPORT ==========\n")
cat(sprintf("1. edges_absr042: n_pairs=%d (expected 751,959)  n_genes=%d (expected 5,450)\n",
            n_pairs, n_genes))
cat(sprintf("2. pickSoftThreshold: power=%d  R²=%.3f  R2>=0.80=%s\n",
            selected_power, r2_selected, if (r2_threshold_met) "YES" else "NO (elbow)"))
if (!is.null(wgcna_labels_out)) {
  cat(sprintf("3. WGCNA: n_modules=%d  grey_rate=%.1f%%  size min/med/max=%d/%.0f/%d  power=%d\n",
              wgcna_n_mods, 100 * wgcna_grey_rate,
              min(sizes_w), median(sizes_w), max(sizes_w), wgcna_power))
} else {
  cat("3. WGCNA: FAILED\n")
}
cat(sprintf("4. Louvain: n_modules=%d  grey_rate=%.1f%%  size min/med/max=%d/%.0f/%d\n",
            lv_n_mods, 100 * lv_grey_rate,
            min(sizes_lv), median(sizes_lv), max(sizes_lv)))
cat(sprintf("5. GO: WGCNA=%s  Louvain=%s modules with ≥1 significant term\n",
            if (is.na(n_go_wgcna)) "SKIP" else as.character(n_go_wgcna),
            if (is.na(n_go_louvain)) "SKIP" else as.character(n_go_louvain)))
if (!is.null(cp_wgcna) && "dominant_pattern" %in% names(cp_wgcna)) {
  dom_w <- sort(table(cp_wgcna$dominant_pattern), decreasing = TRUE)
  cat("6. Condition patterns WGCNA — dominant:", paste(names(dom_w), collapse = " > "), "\n")
}
if ("dominant_pattern" %in% names(cp_louvain)) {
  dom_lv <- sort(table(cp_louvain$dominant_pattern), decreasing = TRUE)
  cat("   Condition patterns Louvain — dominant:", paste(names(dom_lv), collapse = " > "), "\n")
}
cat(sprintf("7. BON3 partners at |r|>=%.2f: %d\n", MIN_ABS_R, nrow(bon3_042)))
cat(sprintf("8. WRKY IDs: %d\n", length(wrky_ids)))
cat(sprintf("9. Total wall time: %.1f min\n", t_total / 60))
cat("=== DONE ===\n")
