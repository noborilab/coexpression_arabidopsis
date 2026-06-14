# inst/scripts/benchmark_modules_pathogen.R
#
# Comprehensive module-construction benchmark — pathogen multiome GGM data.
# Evaluation is STRUCTURE-ONLY: no biological judgment, no curated gene sets,
# no GO-based ranking. GO and TF annotations are computed mechanically for ALL
# modules and saved as reference output only.
#
# Outputs: results/pathogen_multiome/method_benchmark/
#
# See docs/ARCHITECTURE.md and docs/PIPELINE_FLAGS.md for project context.

suppressPackageStartupMessages({
  library(coexpressionArabidopsis)
})

SCRIPT_START <- proc.time()

# ============================================================
# Step 0 — Preliminaries
# ============================================================

RESULTS_DIR <- "results/pathogen_multiome"
BENCH_DIR   <- file.path(RESULTS_DIR, "method_benchmark")
dir.create(BENCH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BENCH_DIR, "assignments"), showWarnings = FALSE)
dir.create(file.path(BENCH_DIR, "profiles"),    showWarnings = FALSE)

message("[", Sys.time(), "] Loading robustness result (not recomputed)...")
rob <- readRDS(file.path(RESULTS_DIR, "robustness", "robustness_result.rds"))
message("  pair_scores rows: ", nrow(rob$pair_scores))

message("[", Sys.time(), "] Loading per-condition network results...")
network_list <- load_network_results(
  "output_per_condition",
  strata = c("Mock", "DC3000", "AvrRpt2", "AvrRpm1"),
  mode   = "singlecellggm"
)
for (s in names(network_list))
  message("  ", s, ": ", nrow(network_list[[s]]$edge_table), " edges")

# igraph
for (pkg in c("igraph")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(igraph)
message("igraph version: ", as.character(packageVersion("igraph")))

# Local null-coalesce (not exported from package)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ============================================================
# Step 1 — Gene symbol map from TAIR10 GFF3
# ============================================================

GFF3_PATH <- paste0("/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/",
                    "SALK_clowd/At_reference/Arabidopsis_thaliana.TAIR10.52.gff3")

message("[", Sys.time(), "] Building gene symbol map from GFF3...")
symbol_map <- tryCatch(
  build_symbol_map(GFF3_PATH),
  error = function(e) {
    message("  WARNING: build_symbol_map() failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(symbol_map)) {
  network_genes <- unique(c(rob$pair_scores$gene_id_A, rob$pair_scores$gene_id_B))
  n_total   <- length(network_genes)
  n_covered <- sum(network_genes %in%
                   symbol_map$gene_id[!is.na(symbol_map$gene_symbol)])
  message("  GFF3 genes: ", nrow(symbol_map))
  message("  Network genes with symbol: ", n_covered, "/", n_total,
          " (", round(100 * n_covered / n_total, 1), "%)")
  write.csv(symbol_map, file.path(RESULTS_DIR, "symbol_map.csv"), row.names = FALSE)
  message("  Saved: symbol_map.csv")
} else {
  message("  Proceeding without symbol map (gene_symbol will be NA everywhere).")
}

lookup_symbol <- function(ids) {
  if (is.null(symbol_map)) return(rep(NA_character_, length(ids)))
  symbol_map$gene_symbol[match(ids, symbol_map$gene_id)]
}

# ============================================================
# Step 2 — Resolve TF metadata path
# ============================================================

TF_CANDIDATES <- c(
  paste0("/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/",
         "TSL/from_Ben/for_tatsuya/data/motifs-2026/Athaliana_motifs_metadata.tsv")
)
hits <- tryCatch(
  list.files(
    paste0("/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI"),
    pattern    = "Athaliana_motifs_metadata",
    recursive  = TRUE,
    full.names = TRUE
  ),
  error = function(e) character(0)
)
message("[", Sys.time(), "] TF metadata search hits: ", length(hits))
if (length(hits) > 0) print(hits)

TF_META_PATH <- NA_character_
for (p in c(TF_CANDIDATES, hits)) {
  if (!is.na(p) && nchar(p) > 0 && file.exists(p)) {
    TF_META_PATH <- p
    break
  }
}
if (is.na(TF_META_PATH)) {
  message("  TF metadata not found — TF annotation will be skipped.")
} else {
  message("  TF metadata resolved: ", TF_META_PATH)
}

# ============================================================
# Step 3 — Benchmark grid
# ============================================================

THRESHOLDS      <- c(0.3, 0.4, 0.5, 0.6, 0.7)
MIN_MODULE_SIZE <- 30L

METHODS <- list(
  list(name = "wgcna_p1", type = "wgcna",  power = 1L),
  list(name = "wgcna_p4", type = "wgcna",  power = 4L),
  list(name = "wgcna_p6", type = "wgcna",  power = 6L),
  list(name = "wgcna_p8", type = "wgcna",  power = 8L),
  list(name = "louvain",  type = "louvain"),
  list(name = "leiden",   type = "leiden")
)
METHOD_NAMES <- vapply(METHODS, `[[`, character(1), "name")
N_METHODS    <- length(METHODS)

message("[", Sys.time(), "] Grid: ", length(THRESHOLDS), " thresholds x ",
        N_METHODS, " methods = ", length(THRESHOLDS) * N_METHODS, " cells")

# Pre-build per-threshold igraph and filtered pair tables
threshold_data <- list()
for (thr in THRESHOLDS) {
  key <- as.character(thr)
  ps  <- rob$pair_scores[!is.na(rob$pair_scores$R_score) &
                         rob$pair_scores$R_score >= thr, , drop = FALSE]
  if (nrow(ps) == 0L) {
    message("  thr=", thr, ": 0 edges — skipped")
    threshold_data[[key]] <- NULL
    next
  }
  gene_ids <- sort(union(ps$gene_id_A, ps$gene_id_B))
  weights  <- pmin(pmax(abs(tanh(ps$z_bar)), 0), 1)

  g <- igraph::graph_from_data_frame(
    d        = data.frame(from = ps$gene_id_A, to = ps$gene_id_B, weight = weights),
    directed = FALSE,
    vertices = data.frame(name = gene_ids)
  )
  threshold_data[[key]] <- list(pairs = ps, gene_ids = gene_ids, graph = g)
  message("  thr=", thr, ": ", nrow(ps), " edges, ", length(gene_ids), " genes")
}

# -----------
# Grid loop
# -----------

assignments <- list()   # assignments[[thr_key]][[method_name]] = data.frame(gene_id, top_module)

for (thr in THRESHOLDS) {
  key <- as.character(thr)
  td  <- threshold_data[[key]]
  if (is.null(td)) next
  assignments[[key]] <- list()

  for (meth in METHODS) {
    label <- paste0("thr", thr, "_", meth$name)
    message("[", Sys.time(), "] ", label)
    out_path <- file.path(BENCH_DIR, "assignments", paste0(label, ".csv"))

    gm <- tryCatch({
      if (meth$type == "wgcna") {
        mi <- build_wgcna_modules(
          rob             = rob,
          network_list    = list(),
          r_score_min     = thr,
          soft_power      = meth$power,
          merge_cut       = 0.25,
          min_module_size = MIN_MODULE_SIZE,
          sub_merge_cut   = NULL
        )
        mi$gene_module[, c("gene_id", "top_module")]

      } else {
        g <- td$graph
        mem <- if (meth$type == "louvain") {
          cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
          igraph::membership(cl)
        } else {
          cl <- igraph::cluster_leiden(
            g,
            weights            = igraph::E(g)$weight,
            objective_function = "modularity"
          )
          igraph::membership(cl)
        }

        gene_ids <- igraph::V(g)$name
        mods     <- as.integer(mem)

        # Map communities smaller than MIN_MODULE_SIZE to module 0
        sz          <- table(mods)
        small_mods  <- as.integer(names(sz)[sz < MIN_MODULE_SIZE])
        mods[mods %in% small_mods] <- 0L

        # Remap non-zero IDs to sequential 1, 2, 3, ...
        nonzero     <- sort(unique(mods[mods > 0L]))
        remap       <- setNames(seq_along(nonzero), nonzero)
        mods[mods > 0L] <- remap[as.character(mods[mods > 0L])]

        data.frame(gene_id = gene_ids, top_module = mods,
                   stringsAsFactors = FALSE)
      }
    }, error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    })

    if (!is.null(gm)) {
      write.csv(gm, out_path, row.names = FALSE)
      message("  assigned=", sum(gm$top_module > 0L),
              " grey=", sum(gm$top_module == 0L),
              " modules=", length(unique(gm$top_module[gm$top_module > 0L])))
    } else {
      message("  FAILED — recording NA metrics")
    }
    assignments[[key]][[meth$name]] <- gm
  }
}

# ============================================================
# Step 4 — Structural metrics (NO biology)
# ============================================================

message("[", Sys.time(), "] Computing structural metrics...")

.gini <- function(x) {
  x <- sort(x[!is.na(x) & x > 0])
  n <- length(x)
  if (n == 0L) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1L) / n
}

metrics_rows <- list()

for (thr in THRESHOLDS) {
  key <- as.character(thr)
  td  <- threshold_data[[key]]
  if (is.null(td)) next

  g       <- td$graph
  ps      <- td$pairs
  edge_df <- igraph::as_data_frame(g, what = "edges")

  for (meth in METHODS) {
    gm <- assignments[[key]][[meth$name]]

    base <- list(
      threshold     = thr,
      method        = meth$name,
      n_input_edges = nrow(ps),
      n_input_genes = length(td$gene_ids)
    )

    if (is.null(gm)) {
      base[c("n_modules","n_assigned","n_grey","pct_grey",
             "module_size_min","module_size_median","module_size_max",
             "module_size_mean","module_size_gini","modularity",
             "mean_intra_weight","mean_inter_weight","separation_ratio")] <- NA
      metrics_rows[[length(metrics_rows) + 1L]] <- as.data.frame(base)
      next
    }

    n_assigned <- sum(gm$top_module > 0L)
    n_grey     <- sum(gm$top_module == 0L)
    n_total    <- nrow(gm)
    mod_sizes  <- as.integer(table(gm$top_module[gm$top_module > 0L]))
    n_modules  <- length(mod_sizes)

    base$n_modules          <- n_modules
    base$n_assigned         <- n_assigned
    base$n_grey             <- n_grey
    base$pct_grey           <- if (n_total > 0L) 100 * n_grey / n_total else NA_real_
    base$module_size_min    <- if (n_modules > 0L) min(mod_sizes)    else NA_real_
    base$module_size_median <- if (n_modules > 0L) median(mod_sizes) else NA_real_
    base$module_size_max    <- if (n_modules > 0L) max(mod_sizes)    else NA_real_
    base$module_size_mean   <- if (n_modules > 0L) mean(mod_sizes)   else NA_real_
    base$module_size_gini   <- if (n_modules > 0L) .gini(mod_sizes)  else NA_real_

    # Modularity computed uniformly on the same input graph for all methods
    base$modularity <- tryCatch({
      vnames  <- igraph::V(g)$name
      mem_vec <- gm$top_module[match(vnames, gm$gene_id)]
      mem_vec[is.na(mem_vec)] <- 0L
      igraph::modularity(g, membership = mem_vec, weights = igraph::E(g)$weight)
    }, error = function(e) NA_real_)

    # Intra / inter mean weight (excluding grey)
    mod_lk <- setNames(gm$top_module, gm$gene_id)
    mod_A  <- mod_lk[edge_df$from]
    mod_B  <- mod_lk[edge_df$to]
    intra  <- !is.na(mod_A) & !is.na(mod_B) & mod_A == mod_B & mod_A > 0L
    inter  <- !is.na(mod_A) & !is.na(mod_B) & mod_A != mod_B & mod_A > 0L & mod_B > 0L

    base$mean_intra_weight <- if (any(intra)) mean(edge_df$weight[intra]) else NA_real_
    base$mean_inter_weight <- if (any(inter)) mean(edge_df$weight[inter]) else NA_real_
    base$separation_ratio  <- if (!is.na(base$mean_inter_weight) &&
                                  base$mean_inter_weight > 0)
      base$mean_intra_weight / base$mean_inter_weight else NA_real_

    metrics_rows[[length(metrics_rows) + 1L]] <- as.data.frame(base)
  }
}

metrics_table <- do.call(rbind, metrics_rows)
write.csv(metrics_table,
          file.path(BENCH_DIR, "structural_metrics.csv"),
          row.names = FALSE)
message("[", Sys.time(), "] Saved: structural_metrics.csv (", nrow(metrics_table), " rows)")

# ============================================================
# Step 5 — Cross-method agreement per threshold
# ============================================================

message("[", Sys.time(), "] Computing cross-method agreement...")

.ari <- function(x, y) {
  ct  <- table(x, y)
  n   <- sum(ct)
  if (n == 0L) return(NA_real_)
  ch2 <- function(v) sum(v * (v - 1L)) / 2
  num <- ch2(as.vector(ct)) - ch2(rowSums(ct)) * ch2(colSums(ct)) / ch2(n)
  den <- (ch2(rowSums(ct)) + ch2(colSums(ct))) / 2 -
         ch2(rowSums(ct)) * ch2(colSums(ct)) / ch2(n)
  if (den == 0) return(0)
  num / den
}

.nmi <- function(x, y) {
  ct   <- table(x, y)
  n    <- sum(ct)
  if (n == 0L) return(NA_real_)
  p_ij <- ct / n
  p_i  <- rowSums(p_ij)
  p_j  <- colSums(p_ij)
  H_x  <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
  H_y  <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))
  denom <- (H_x + H_y) / 2
  if (denom <= 0) return(0)
  MI <- sum(ifelse(p_ij > 0 & outer(p_i, p_j) > 0,
                   p_ij * log(p_ij / outer(p_i, p_j)), 0))
  MI / denom
}

agree_rows <- list()

for (thr in THRESHOLDS) {
  key <- as.character(thr)
  if (is.null(threshold_data[[key]])) next

  assign_at_thr <- assignments[[key]]

  ari_mat <- matrix(NA_real_, nrow = N_METHODS, ncol = N_METHODS,
                    dimnames = list(METHOD_NAMES, METHOD_NAMES))
  diag(ari_mat) <- 1.0

  idx_pairs <- combn(N_METHODS, 2)
  for (k in seq_len(ncol(idx_pairs))) {
    i  <- idx_pairs[1, k]
    j  <- idx_pairs[2, k]
    ni <- METHOD_NAMES[i]
    nj <- METHOD_NAMES[j]
    gi <- assign_at_thr[[ni]]
    gj <- assign_at_thr[[nj]]

    if (is.null(gi) || is.null(gj)) {
      agree_rows[[length(agree_rows) + 1L]] <- data.frame(
        threshold = thr, method_1 = ni, method_2 = nj,
        n_common_assigned = NA_integer_, ARI = NA_real_, NMI = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    common <- intersect(gi$gene_id[gi$top_module > 0L],
                        gj$gene_id[gj$top_module > 0L])
    n_common <- length(common)

    if (n_common < 2L) {
      ari_val <- NA_real_; nmi_val <- NA_real_
    } else {
      x <- gi$top_module[match(common, gi$gene_id)]
      y <- gj$top_module[match(common, gj$gene_id)]
      ari_val <- tryCatch(.ari(x, y), error = function(e) NA_real_)
      nmi_val <- tryCatch(.nmi(x, y), error = function(e) NA_real_)
    }

    ari_mat[i, j] <- ari_val
    ari_mat[j, i] <- ari_val

    agree_rows[[length(agree_rows) + 1L]] <- data.frame(
      threshold         = thr,
      method_1          = ni,
      method_2          = nj,
      n_common_assigned = n_common,
      ARI               = ari_val,
      NMI               = nmi_val,
      stringsAsFactors  = FALSE
    )
  }

  # Per-threshold ARI matrix CSV
  ari_df <- as.data.frame(ari_mat)
  ari_df <- cbind(method = rownames(ari_df), ari_df)
  write.csv(ari_df,
            file.path(BENCH_DIR, paste0("ari_matrix_thr", thr, ".csv")),
            row.names = FALSE)
  message("  thr=", thr, ": ARI matrix saved")
}

if (length(agree_rows) > 0L) {
  agree_tbl <- do.call(rbind, agree_rows)
  write.csv(agree_tbl, file.path(BENCH_DIR, "cross_method_agreement.csv"),
            row.names = FALSE)
  message("[", Sys.time(), "] Saved: cross_method_agreement.csv")
}

# ============================================================
# Step 6 — Per-module reference profiles at threshold = 0.5
# ============================================================

PROFILE_THR <- 0.5
prof_key    <- as.character(PROFILE_THR)

message("[", Sys.time(), "] Building reference profiles at threshold = ", PROFILE_THR, "...")

# GO: all significant terms saved as a long table (reference output only)
.go_long <- function(gene_mod, universe, org_db = "org.At.tair.db", pval_cut = 0.05) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace(org_db, quietly = TRUE)) {
    message("  clusterProfiler or ", org_db, " not available; skipping GO.")
    return(NULL)
  }
  db   <- get(org_db, envir = asNamespace(org_db))
  mods <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))
  rows <- list()
  for (m in mods) {
    tryCatch({
      ego <- clusterProfiler::enrichGO(
        gene          = gene_mod$gene_id[gene_mod$top_module == m],
        universe      = universe,
        OrgDb         = db,
        keyType       = "TAIR",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cut,
        readable      = FALSE
      )
      if (!is.null(ego)) {
        res <- ego@result[!is.na(ego@result$p.adjust) &
                          ego@result$p.adjust < pval_cut, , drop = FALSE]
        if (nrow(res) > 0L)
          rows[[length(rows) + 1L]] <- data.frame(
            module_id  = as.integer(m),
            GO_id      = res$ID,
            term       = res$Description,
            p.adjust   = res$p.adjust,
            gene_count = res$Count,
            stringsAsFactors = FALSE
          )
      }
    }, error = function(e)
      message("  GO error module ", m, ": ", conditionMessage(e)))
  }
  if (length(rows) > 0L) do.call(rbind, rows) else NULL
}

# Condition profiles: per-module per-condition mean intra-module edge weight
.cond_profiles <- function(gene_mod, network_list) {
  conds <- names(network_list)
  mods  <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))
  rows  <- lapply(mods, function(m) {
    mg <- gene_mod$gene_id[gene_mod$top_module == m]
    lapply(conds, function(cond) {
      et    <- network_list[[cond]]$edge_table
      intra <- et[et$gene_id_A %in% mg & et$gene_id_B %in% mg, , drop = FALSE]
      data.frame(module_id = as.integer(m), condition = cond,
                 mean_intra_weight = if (nrow(intra) > 0L) mean(abs(intra$weight)) else 0,
                 stringsAsFactors = FALSE)
    })
  })
  result <- do.call(rbind, lapply(rows, function(x) do.call(rbind, x)))
  if ("Mock" %in% result$condition) {
    mock_w <- setNames(
      result$mean_intra_weight[result$condition == "Mock"],
      as.character(result$module_id[result$condition == "Mock"])
    )
    result$delta_vs_mock <- result$mean_intra_weight -
      mock_w[as.character(result$module_id)]
  } else {
    result$delta_vs_mock <- NA_real_
  }
  result
}

# Hub genes: top 20 per module by intra-module node strength
.hub_genes <- function(gene_mod, g, n_top = 20L) {
  edge_df <- igraph::as_data_frame(g, what = "edges")
  mods    <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))
  rows    <- lapply(mods, function(m) {
    mg    <- gene_mod$gene_id[gene_mod$top_module == m]
    intra <- edge_df[edge_df$from %in% mg & edge_df$to %in% mg, , drop = FALSE]
    if (nrow(intra) == 0L) return(NULL)
    from_s <- tapply(intra$weight, intra$from, sum)
    to_s   <- tapply(intra$weight, intra$to,   sum)
    all_g  <- union(names(from_s), names(to_s))
    str_v  <- vapply(all_g, function(id)
      (from_s[id] %||% 0) + (to_s[id] %||% 0), numeric(1))
    str_v  <- sort(str_v[names(str_v) %in% mg], decreasing = TRUE)
    top_g  <- head(names(str_v), n_top)
    data.frame(module_id   = as.integer(m),
               gene_id     = top_g,
               gene_symbol = lookup_symbol(top_g),
               hub_rank    = seq_along(top_g),
               hub_score   = str_v[top_g],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

if (!is.null(threshold_data[[prof_key]])) {
  universe_genes <- threshold_data[[prof_key]]$gene_ids
  g_prof         <- threshold_data[[prof_key]]$graph

  for (meth in METHODS) {
    mname <- meth$name
    gm    <- assignments[[prof_key]][[mname]]

    if (is.null(gm)) {
      message("  Skipping profiles for ", mname, " (no assignment)")
      next
    }

    pdir <- file.path(BENCH_DIR, "profiles", mname)
    dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

    n_mods <- length(unique(gm$top_module[gm$top_module > 0L]))
    message("[", Sys.time(), "] Profile: ", mname, " — ", n_mods, " modules")

    # gene_module.csv
    write.csv(gm, file.path(pdir, "gene_module.csv"), row.names = FALSE)

    # module_meta skeleton
    mod_counts <- table(gm$top_module[gm$top_module > 0L])
    module_meta <- data.frame(
      module_id = as.integer(names(mod_counts)),
      n_genes   = as.integer(mod_counts),
      stringsAsFactors = FALSE
    )

    # Condition profiles
    tryCatch({
      cp <- .cond_profiles(gm, network_list)
      write.csv(cp, file.path(pdir, "module_condition_profiles.csv"), row.names = FALSE)
      # Top condition per module
      for (m in module_meta$module_id) {
        sub <- cp[cp$module_id == m, ]
        module_meta$top_condition[module_meta$module_id == m] <-
          sub$condition[which.max(sub$mean_intra_weight)]
      }
    }, error = function(e)
      message("  Condition profiles error: ", conditionMessage(e)))

    # Hub genes
    tryCatch({
      hg <- .hub_genes(gm, g_prof)
      if (!is.null(hg) && nrow(hg) > 0L)
        write.csv(hg, file.path(pdir, "hub_genes.csv"), row.names = FALSE)
    }, error = function(e)
      message("  Hub genes error: ", conditionMessage(e)))

    # TF annotation
    if (!is.na(TF_META_PATH)) {
      tryCatch({
        mi_tmp <- list(
          gene_module = gm,
          module_meta = module_meta,
          module_tfs  = data.frame(module_id = integer(), gene_id = character(),
                                   gene_symbol = character(), tf_family = character(),
                                   stringsAsFactors = FALSE)
        )
        mi_tmp <- annotate_tfs(mi_tmp, TF_META_PATH)
        if (!is.null(mi_tmp$module_tfs) && nrow(mi_tmp$module_tfs) > 0L)
          mi_tmp$module_tfs$gene_symbol <- lookup_symbol(mi_tmp$module_tfs$gene_id)
        write.csv(mi_tmp$module_tfs, file.path(pdir, "module_tfs.csv"), row.names = FALSE)
      }, error = function(e)
        message("  TF annotation error: ", conditionMessage(e)))
    } else {
      write.csv(
        data.frame(module_id=integer(), gene_id=character(),
                   gene_symbol=character(), tf_family=character(),
                   stringsAsFactors=FALSE),
        file.path(pdir, "module_tfs.csv"), row.names = FALSE
      )
    }

    # GO enrichment — all significant terms, saved as long table (reference only)
    tryCatch({
      go_long <- .go_long(gm, universe_genes)
      if (!is.null(go_long) && nrow(go_long) > 0L) {
        write.csv(go_long, file.path(pdir, "go_enrichment_long.csv"), row.names = FALSE)
        message("  GO: ", nrow(go_long), " significant term rows")
        # Populate go_top in module_meta (top term by p.adjust, reference only)
        for (m in module_meta$module_id) {
          sub <- go_long[go_long$module_id == m, ]
          if (nrow(sub) > 0L)
            module_meta$go_top[module_meta$module_id == m] <-
              sub$term[which.min(sub$p.adjust)]
        }
      } else {
        write.csv(
          data.frame(module_id=integer(), GO_id=character(), term=character(),
                     p.adjust=numeric(), gene_count=integer(), stringsAsFactors=FALSE),
          file.path(pdir, "go_enrichment_long.csv"), row.names = FALSE
        )
        message("  GO: no significant terms (or clusterProfiler unavailable)")
      }
    }, error = function(e)
      message("  GO error: ", conditionMessage(e)))

    # Preservation fallback (reference only, not used for ranking)
    tryCatch({
      mi_tmp <- list(gene_module = gm)
      pres   <- compute_preservation_fallback(mi_tmp, network_list, n_perm = 100L)
      if (!is.null(pres) && nrow(pres) > 0L) {
        m_idx <- match(module_meta$module_id, pres$module_id)
        module_meta$zsummary[!is.na(m_idx)] <- pres$zsummary[m_idx[!is.na(m_idx)]]
        module_meta$preservation_method[!is.na(m_idx)] <- "fallback_meancor"
      }
    }, error = function(e)
      message("  Preservation error: ", conditionMessage(e)))

    write.csv(module_meta, file.path(pdir, "module_meta.csv"), row.names = FALSE)
    message("  Done: ", pdir)
  }
} else {
  message("  Threshold 0.5 data unavailable; skipping profiles.")
}

# ============================================================
# Step 7 — Benchmark report (structural facts only)
# ============================================================

message("[", Sys.time(), "] Writing BENCHMARK_REPORT.md...")

# Helper: format data.frame as markdown table
.md_table <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) format(round(x, 3), nsmall = 0) else as.character(x)
  })
  df[is.na(df)] <- "NA"
  hdr  <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep  <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, rows), collapse = "\n")
}

report_lines <- c(
  "# Module-Construction Benchmark Report — Pathogen Multiome",
  "",
  paste0("Run completed: ", Sys.time()),
  "",
  "---",
  "",
  "## 1. Grid Description",
  "",
  paste0("- Thresholds (R_score): ",
         paste(THRESHOLDS, collapse = ", ")),
  paste0("- Methods: ", paste(METHOD_NAMES, collapse = ", ")),
  paste0("- Total cells: ", length(THRESHOLDS) * N_METHODS),
  "",
  paste0("- Cells run: ",
         if (exists("metrics_table"))
           sum(!is.na(metrics_table$n_modules))
         else "—"),
  paste0("- Cells failed (NA): ",
         if (exists("metrics_table"))
           sum(is.na(metrics_table$n_modules))
         else "—"),
  "",
  "---",
  "",
  "## 2. Structural Metrics (all 30 cells)",
  "",
  .md_table(metrics_table),
  "",
  "---",
  "",
  "## 3. Per-Threshold Cross-Method ARI Matrices"
)

for (thr in THRESHOLDS) {
  key  <- as.character(thr)
  path <- file.path(BENCH_DIR, paste0("ari_matrix_thr", thr, ".csv"))
  if (!file.exists(path)) next
  mat_df <- read.csv(path, stringsAsFactors = FALSE)
  report_lines <- c(report_lines,
    "", paste0("### Threshold = ", thr),
    "", .md_table(mat_df))
}

# Structural observations
obs_lines <- c(
  "",
  "---",
  "",
  "## 4. Structural Observations",
  ""
)

if (exists("metrics_table") && nrow(metrics_table) > 0L) {
  mt <- metrics_table[!is.na(metrics_table$modularity), ]

  for (thr in THRESHOLDS) {
    sub <- mt[mt$threshold == thr, ]
    if (nrow(sub) == 0L) next

    max_mod_m <- sub$method[which.max(sub$modularity)]
    min_grey_m <- sub$method[which.min(sub$pct_grey)]
    obs_lines <- c(obs_lines,
      paste0("**Threshold ", thr, ":**"),
      paste0("- n_input_edges=", sub$n_input_edges[1],
             ", n_input_genes=", sub$n_input_genes[1]),
      paste0("- Modularity range: ",
             round(min(sub$modularity, na.rm = TRUE), 3), " – ",
             round(max(sub$modularity, na.rm = TRUE), 3),
             " (highest: ", max_mod_m, ")"),
      paste0("- Grey rate range: ",
             round(min(sub$pct_grey, na.rm = TRUE), 1), "% – ",
             round(max(sub$pct_grey, na.rm = TRUE), 1), "%",
             " (lowest grey: ", min_grey_m, ")"),
      paste0("- Module count range: ",
             min(sub$n_modules, na.rm = TRUE), " – ",
             max(sub$n_modules, na.rm = TRUE)),
      paste0("- Median module size range: ",
             round(min(sub$module_size_median, na.rm = TRUE), 0), " – ",
             round(max(sub$module_size_median, na.rm = TRUE), 0)),
      ""
    )
  }

  if (exists("agree_tbl")) {
    for (thr in THRESHOLDS) {
      sub <- agree_tbl[agree_tbl$threshold == thr & !is.na(agree_tbl$ARI), ]
      if (nrow(sub) == 0L) next
      max_pair <- sub[which.max(sub$ARI), ]
      min_pair <- sub[which.min(sub$ARI), ]
      obs_lines <- c(obs_lines,
        paste0("**Cross-method ARI at threshold ", thr, ":**"),
        paste0("- Highest pair: ", max_pair$method_1, " vs ", max_pair$method_2,
               " (ARI=", round(max_pair$ARI, 3), ", n_common=", max_pair$n_common_assigned, ")"),
        paste0("- Lowest pair: ", min_pair$method_1, " vs ", min_pair$method_2,
               " (ARI=", round(min_pair$ARI, 3), ", n_common=", min_pair$n_common_assigned, ")"),
        ""
      )
    }
  }
}

design_lines <- c(
  "---",
  "",
  "## 5. Design Note: benchmark_module_methods() generalization",
  "",
  "This benchmark could become a reusable `benchmark_module_methods()` function",
  "within the pipeline. Proposed inputs/outputs:",
  "",
  "**Inputs:**",
  "- `rob`: RobustnessResult (or any named edge list with weight column)",
  "- `thresholds`: numeric vector of filter thresholds",
  "- `methods`: a method registry (list of named specs: type, params)",
  "- `min_module_size`: integer",
  "- `network_list`: for annotation steps (optional)",
  "",
  "**Outputs:**",
  "- `structural_metrics.csv` (all cells: counts, modularity, grey rate, separation)",
  "- `cross_method_agreement.csv` (pairwise ARI/NMI per threshold)",
  "- `ari_matrix_thr{X}.csv` (per-threshold matrices)",
  "- `assignments/{cell}.csv` (gene-level assignments, all cells)",
  "- `profiles/{method}/` (annotated module tables at a chosen threshold)",
  "- `BENCHMARK_REPORT.md`",
  "",
  "**What needs generalizing:**",
  "- Method registry: currently hard-coded; should be a config-driven list",
  "  with type dispatch (wgcna/louvain/leiden + future methods)",
  "- Modularity computation: already uniform across methods (igraph) — keep",
  "- Annotation steps (GO, TF, preservation) should be optional flags, since",
  "  they are slow and irrelevant for pure structural evaluation",
  "- Profile threshold should be a configurable parameter, not hard-coded to 0.5",
  "- The function should return a structured list of all metric tables for",
  "  programmatic access, not just write files",
  "",
  "**Observation:** GGM output is a sparse partial-correlation graph. WGCNA",
  "soft-thresholding was designed for dense correlation matrices; auto power",
  "selection on GGM input chose power=1 and scale-free fit was poor (max R²≈0.67).",
  "Graph-clustering methods (Louvain/Leiden) operate directly on weighted graphs",
  "without a soft-thresholding step, which may be structurally more natural for",
  "GGM output. The benchmark quantifies this difference via modularity and grey rate",
  "without asserting which is biologically preferable."
)

writeLines(
  c(report_lines, obs_lines, design_lines),
  file.path(BENCH_DIR, "BENCHMARK_REPORT.md")
)
message("[", Sys.time(), "] Saved: BENCHMARK_REPORT.md")

# ============================================================
# Step 8 — Update docs (appended in-place)
# ============================================================

flag11 <- paste0(
  "\n---\n\n",
  "## FLAG-11: WGCNA soft-thresholding on GGM partial-correlation networks\n",
  "**Phase**: benchmark (2026-06)\n",
  "**Issue**: WGCNA soft-thresholding was designed for dense correlation matrices.\n",
  "On GGM output (already sparse; partial correlations), auto soft-power selection\n",
  "chose power=1 and scale-free R² fit was poor (max R²≈0.67). The benchmark\n",
  "(inst/scripts/benchmark_modules_pathogen.R) compares WGCNA at powers 1/4/6/8\n",
  "and graph-clustering methods (Louvain, Leiden) across five R_score thresholds\n",
  "(0.3–0.7) on the pathogen multiome GGM robustness results. Evaluation is\n",
  "structure-only: modularity, grey rate, module size, cross-method ARI.\n",
  "See results/pathogen_multiome/method_benchmark/BENCHMARK_REPORT.md for findings.\n",
  "**Decision**: Under evaluation. User to decide method based on benchmark\n",
  "structural summary plus biological review.\n",
  "**Status**: Benchmark running; results pending.\n"
)

arch_section <- paste0(
  "\n---\n\n",
  "## Module construction methods (under evaluation)\n\n",
  "For GGM mode, the interpretation layer currently supports multiple module\n",
  "construction paths:\n\n",
  "- **WGCNA** at explicit soft powers (1, 4, 6, 8) — `build_wgcna_modules()`\n",
  "  with `soft_power` override. Auto power selection on GGM output converges\n",
  "  to power=1 with poor scale-free fit (FLAG-11).\n",
  "- **Louvain** — `igraph::cluster_louvain()` on the R_score-filtered weighted\n",
  "  graph (weight = abs(tanh(z_bar))).\n",
  "- **Leiden** — `igraph::cluster_leiden()` with objective_function='modularity'\n",
  "  on the same graph.\n\n",
  "A comprehensive benchmarking layer (5 thresholds × 6 methods = 30 cells)\n",
  "evaluates structure-only metrics (modularity, grey rate, module sizes, Gini,\n",
  "cross-method ARI/NMI) without biological pre-judgment.\n",
  "See `inst/scripts/benchmark_modules_pathogen.R` and\n",
  "`results/pathogen_multiome/method_benchmark/`.\n\n",
  "This benchmarking step may become a permanent, general-purpose\n",
  "`benchmark_module_methods()` pipeline component — see BENCHMARK_REPORT.md\n",
  "design note section for the proposed API.\n"
)

tryCatch({
  cat(flag11, file = "docs/PIPELINE_FLAGS.md", append = TRUE)
  message("[", Sys.time(), "] Appended FLAG-11 to docs/PIPELINE_FLAGS.md")
}, error = function(e)
  message("  Doc update failed: ", conditionMessage(e)))

tryCatch({
  cat(arch_section, file = "docs/ARCHITECTURE.md", append = TRUE)
  message("[", Sys.time(), "] Appended module-construction section to docs/ARCHITECTURE.md")
}, error = function(e)
  message("  Doc update failed: ", conditionMessage(e)))

# ============================================================
# Done
# ============================================================

elapsed <- (proc.time() - SCRIPT_START)[["elapsed"]]
message("[", Sys.time(), "] COMPLETE. Wall time: ",
        sprintf("%.1f min (%.0f sec)", elapsed / 60, elapsed))
