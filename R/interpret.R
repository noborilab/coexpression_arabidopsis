#' @title Module Construction and Biological Interpretation
#'
#' @description
#' Shared output/interpretation layer used by **both** estimation modes.
#'
#' Covers:
#' - WGCNA-style signed-network module construction (soft power, merge threshold)
#' - Granularity hierarchy: coarse top-level modules with nested sub-modules
#' - Cross-context module preservation (lightweight mean intramodular |cor|
#'   z-score fallback — `modulePreservation` timed out on the pathogen side
#'   in the CZL run on ~1,500 samples × 11k genes; FLAG-02)
#' - Hub genes (kME)
#' - TF intersection (regulator hints per module; use lab TF file, NOT PlantTFDB)
#' - GO BP enrichment (clusterProfiler, graceful skip if not installed)
#'
#' @name interpret
NULL

# ---------------------------------------------------------------------------
# build_wgcna_modules
# ---------------------------------------------------------------------------

#' Build WGCNA co-expression modules from a robustness-filtered edge set
#'
#' Constructs an adjacency matrix from the robustness-weighted mean partial
#' correlations, applies WGCNA soft-thresholding, TOM, and hierarchical
#' clustering to produce coarse top-level modules and optional finer
#' sub-modules.
#'
#' @param rob RobustnessResult from [compute_robustness()].
#' @param network_list Named list of NetworkResult (same strata as `rob`).
#'   Used by downstream annotation functions (context, preservation).
#' @param r_score_min Minimum R_score to include an edge. Default `0` (all
#'   edges). Ignored when `min_abs_r` is supplied.
#' @param min_abs_r When non-`NULL`, filter edges by mean |Spearman r| =
#'   `abs(tanh(z_bar)) >= min_abs_r` instead of R_score. This is the
#'   preferred sparsification lever for pseudobulk networks because R_score
#'   is a discrete condition-count (1–4) that leaves the network very dense
#'   even at R_score = 1.0. Default `NULL` (use `r_score_min`).
#' @param soft_power WGCNA soft-thresholding power. `NULL` (default) =
#'   auto-pick by scale-free fit with `pickSoftThreshold.fromSimilarity`.
#' @param merge_cut Module merge threshold (height). Default `0.25`. Lower
#'   value = fewer, coarser modules.
#' @param min_module_size Minimum number of genes per module. Default `30`.
#' @param sub_merge_cut Merge threshold for sub-modules (finer granularity).
#'   Default `0.10`. `NULL` = skip sub-module construction.
#' @return ModuleInput (named list; see `docs/OUTPUT_SCHEMA.md`).
#' @export
build_wgcna_modules <- function(rob,
                                network_list,
                                r_score_min     = 0,
                                min_abs_r       = NULL,
                                soft_power      = NULL,
                                merge_cut       = 0.25,
                                min_module_size = 30,
                                sub_merge_cut   = 0.10) {

  if (!requireNamespace("WGCNA", quietly = TRUE))
    stop("WGCNA must be installed. Run: install.packages('WGCNA') or ",
         "BiocManager::install('WGCNA')")

  # 1. Filter pairs -------------------------------------------------------------
  ps <- rob$pair_scores
  if (!is.null(min_abs_r)) {
    ps <- ps[!is.na(ps$z_bar) & abs(tanh(ps$z_bar)) >= min_abs_r, , drop = FALSE]
    if (nrow(ps) == 0L)
      stop("No pairs remain after min_abs_r filter (min_abs_r = ", min_abs_r, ").")
  } else {
    ps <- ps[!is.na(ps$R_score) & ps$R_score >= r_score_min, , drop = FALSE]
    if (nrow(ps) == 0L)
      stop("No pairs remain after r_score_min filter (r_score_min = ", r_score_min, ").")
  }

  # 2. Gene universe ------------------------------------------------------------
  gene_ids <- sort(union(ps$gene_id_A, ps$gene_id_B))
  n_genes  <- length(gene_ids)
  gene_idx <- setNames(seq_along(gene_ids), gene_ids)

  # 3. Adjacency matrix ---------------------------------------------------------
  # Weight = abs(tanh(z_bar)) = back-transformed random-effects mean correlation.
  weights <- abs(tanh(ps$z_bar))
  weights <- pmin(weights, 1.0)

  A <- matrix(0.0, nrow = n_genes, ncol = n_genes,
              dimnames = list(gene_ids, gene_ids))
  i_A <- gene_idx[ps$gene_id_A]
  i_B <- gene_idx[ps$gene_id_B]
  A[cbind(i_A, i_B)] <- weights
  A[cbind(i_B, i_A)] <- weights
  diag(A) <- 1.0

  # 4. Soft-thresholding --------------------------------------------------------
  if (is.null(soft_power)) {
    sft <- tryCatch(
      WGCNA::pickSoftThreshold.fromSimilarity(
        A, powerVector = 1:20, verbose = 0
      ),
      error = function(e) {
        message("pickSoftThreshold.fromSimilarity failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(sft)) {
      fit_tbl <- sft$fitIndices
      good    <- fit_tbl[!is.na(fit_tbl$SFT.R.sq) & fit_tbl$SFT.R.sq >= 0.85, ]
      if (nrow(good) > 0L) {
        soft_power <- good$Power[1L]
      } else {
        soft_power <- fit_tbl$Power[which.max(fit_tbl$SFT.R.sq)]
        message("No power achieved scale-free R² ≥ 0.85; using best available ",
                "(power = ", soft_power, ", R² = ",
                round(max(fit_tbl$SFT.R.sq, na.rm = TRUE), 3), ")")
      }
    } else {
      soft_power <- 6L
      message("Falling back to soft_power = 6")
    }
    message("Auto-selected soft_power = ", soft_power)
  }

  A_soft <- A ^ soft_power

  # 5. TOM and hierarchical clustering ------------------------------------------
  TOM      <- WGCNA::TOMsimilarity(A_soft, TOMType = "signed", verbose = 0)
  dimnames(TOM) <- dimnames(A_soft)
  dist_mat <- 1 - TOM
  tree     <- hclust(as.dist(dist_mat), method = "average")

  # 6. Initial module cut -------------------------------------------------------
  if (!requireNamespace("dynamicTreeCut", quietly = TRUE))
    stop("dynamicTreeCut must be installed: install.packages('dynamicTreeCut')")

  init_labels <- dynamicTreeCut::cutreeDynamic(
    dendro         = tree,
    distM          = dist_mat,
    deepSplit      = 2,
    minClusterSize = min_module_size,
    method         = "hybrid",
    verbose        = 0
  )

  # 7a. Top-level modules (coarse merge) ----------------------------------------
  top_merged <- WGCNA::mergeCloseModules(
    exprData  = t(A_soft),
    colors    = init_labels,
    cutHeight = merge_cut,
    verbose   = 0
  )
  top_labels <- as.integer(top_merged$colors)

  # 7b. Sub-modules (finer merge from same init_labels) -------------------------
  sub_labels <- rep(0L, n_genes)
  if (!is.null(sub_merge_cut)) {
    sub_merged <- WGCNA::mergeCloseModules(
      exprData  = t(A_soft),
      colors    = init_labels,
      cutHeight = sub_merge_cut,
      verbose   = 0
    )
    sub_labels <- as.integer(sub_merged$colors)
  }

  # 8. Module eigengenes --------------------------------------------------------
  MEs <- WGCNA::moduleEigengenes(
    expr = t(A_soft), colors = top_labels, verbose = 0
  )$eigengenes

  # 9. kME = correlation of each gene's adjacency profile with its ME -----------
  kme_mat <- tryCatch(
    cor(t(A_soft), MEs),
    error = function(e) NULL
  )

  kme_vec <- vapply(seq_along(gene_ids), function(i) {
    mod <- top_labels[i]
    if (mod == 0L || is.null(kme_mat)) return(NA_real_)
    me_col <- paste0("ME", mod)
    if (me_col %in% colnames(kme_mat)) kme_mat[gene_ids[i], me_col] else NA_real_
  }, numeric(1))

  # 10. Assemble ModuleInput ----------------------------------------------------

  gene_module <- data.frame(
    gene_id    = gene_ids,
    top_module = top_labels,
    sub_module = sub_labels,
    kME        = kme_vec,
    stringsAsFactors = FALSE
  )

  # module_meta: one row per non-zero module
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

  # module_hier: sub_module → top_module mapping derived from gene assignments
  hier_df <- unique(data.frame(
    sub_module = gene_module$sub_module,
    top_module = gene_module$top_module,
    stringsAsFactors = FALSE
  ))
  module_hier <- hier_df[hier_df$sub_module > 0L, , drop = FALSE]
  module_hier <- module_hier[order(module_hier$sub_module), ]
  rownames(module_hier) <- NULL

  # hub_genes: top 20 per module by kME
  unique_mods <- as.integer(names(mod_counts))
  hub_list <- lapply(unique_mods, function(m) {
    sub_gm <- gene_module[gene_module$top_module == m & !is.na(gene_module$kME), ]
    sub_gm <- sub_gm[order(sub_gm$kME, decreasing = TRUE), ]
    sub_gm <- head(sub_gm, 20L)
    if (nrow(sub_gm) == 0L) return(NULL)
    data.frame(
      module_id   = as.integer(m),
      gene_id     = sub_gm$gene_id,
      gene_symbol = NA_character_,
      kME         = sub_gm$kME,
      hub_rank    = seq_len(nrow(sub_gm)),
      stringsAsFactors = FALSE
    )
  })
  hub_genes <- do.call(rbind, Filter(Negate(is.null), hub_list))
  if (is.null(hub_genes) || nrow(hub_genes) == 0L) {
    hub_genes <- data.frame(
      module_id = integer(), gene_id = character(), gene_symbol = character(),
      kME = numeric(), hub_rank = integer(), stringsAsFactors = FALSE
    )
  }

  # module_tfs: empty — filled by annotate_tfs()
  module_tfs <- data.frame(
    module_id   = integer(),
    gene_id     = character(),
    gene_symbol = character(),
    tf_family   = character(),
    stringsAsFactors = FALSE
  )

  list(
    gene_module  = gene_module,
    module_meta  = module_meta,
    module_hier  = module_hier,
    hub_genes    = hub_genes,
    module_tfs   = module_tfs,
    eigengenes   = MEs
  )
}

# ---------------------------------------------------------------------------
# compute_preservation_fallback  (FLAG-02)
# ---------------------------------------------------------------------------

#' Lightweight module preservation fallback
#'
#' Computes a z-score proxy for module preservation by comparing mean
#' intramodular |weight| in a test network against a permutation null.
#' Used in place of `WGCNA::modulePreservation`, which timed out on the
#' pathogen data (~1,500 samples × 11k genes) in the CZL run (FLAG-02).
#'
#' @param mod_input ModuleInput from [build_wgcna_modules()].
#' @param network_list2 Named list of NetworkResult used as the preservation
#'   test network (a different context from the one used to build modules).
#' @param n_perm Number of permutation draws for the null. Default `100`.
#' @return `data.frame` with columns `module_id`, `zsummary`,
#'   `preservation_method` (always `"fallback_meancor"`).
#' @export
compute_preservation_fallback <- function(mod_input, network_list2,
                                          n_perm = 100L) {

  gene_mod <- mod_input$gene_module

  # Pool all edges from network_list2
  all_edges <- do.call(rbind, lapply(network_list2, function(nr) nr$edge_table))
  if (is.null(all_edges) || nrow(all_edges) == 0L)
    stop("network_list2 contains no edges.")

  all_gene_pool <- unique(c(all_edges$gene_id_A, all_edges$gene_id_B))
  unique_mods   <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))

  rows <- lapply(unique_mods, function(m) {
    mod_genes <- gene_mod$gene_id[gene_mod$top_module == m]
    intra     <- all_edges[all_edges$gene_id_A %in% mod_genes &
                           all_edges$gene_id_B %in% mod_genes, , drop = FALSE]
    obs_mean  <- if (nrow(intra) > 0L) mean(abs(intra$weight)) else 0.0

    n_mod <- length(mod_genes)
    pool_n <- min(n_mod, length(all_gene_pool))
    perm_means <- vapply(seq_len(n_perm), function(.) {
      perm_g <- sample(all_gene_pool, pool_n)
      perm_intra <- all_edges[all_edges$gene_id_A %in% perm_g &
                              all_edges$gene_id_B %in% perm_g, , drop = FALSE]
      if (nrow(perm_intra) > 0L) mean(abs(perm_intra$weight)) else 0.0
    }, numeric(1))

    perm_mean <- mean(perm_means)
    perm_sd   <- sd(perm_means)
    zsummary  <- if (!is.na(perm_sd) && perm_sd > 0) {
      (obs_mean - perm_mean) / perm_sd
    } else 0.0

    data.frame(module_id = as.integer(m), zsummary = zsummary,
               preservation_method = "fallback_meancor",
               stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# annotate_context
# ---------------------------------------------------------------------------

#' Add condition-level context to module_meta
#'
#' For each module, computes the mean intramodular edge weight per condition
#' and identifies the condition with the highest activity. Also computes
#' `delta_treatment` relative to Mock.
#'
#' @param mod_input ModuleInput from [build_wgcna_modules()].
#' @param network_list Named list of NetworkResult (one per condition).
#' @return `mod_input` with `$module_meta$top_organ_or_condition` and
#'   `$module_meta$delta_treatment` filled.
#' @export
annotate_context <- function(mod_input, network_list) {
  gene_mod    <- mod_input$gene_module
  conditions  <- names(network_list)
  unique_mods <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))

  for (m in unique_mods) {
    mod_genes <- gene_mod$gene_id[gene_mod$top_module == m]

    cond_weights <- vapply(conditions, function(cond) {
      et    <- network_list[[cond]]$edge_table
      intra <- et[et$gene_id_A %in% mod_genes & et$gene_id_B %in% mod_genes, ]
      if (nrow(intra) > 0L) mean(abs(intra$weight)) else 0.0
    }, numeric(1))

    top_cond <- conditions[which.max(cond_weights)]

    mock_w <- cond_weights["Mock"]
    delta_str <- if (!is.na(mock_w)) {
      non_mock <- conditions[conditions != "Mock"]
      parts    <- vapply(non_mock, function(cond)
        sprintf("%s:%+.3f", cond, cond_weights[cond] - mock_w),
        character(1))
      paste(parts, collapse = ";")
    } else {
      NA_character_
    }

    ri <- which(mod_input$module_meta$module_id == m)
    mod_input$module_meta$top_organ_or_condition[ri] <- top_cond
    mod_input$module_meta$delta_treatment[ri]        <- delta_str
  }

  mod_input
}

# ---------------------------------------------------------------------------
# annotate_go
# ---------------------------------------------------------------------------

#' GO BP enrichment per module
#'
#' Runs `clusterProfiler::enrichGO()` for each module and fills
#' `module_meta$go_top` with the top term (lowest adjusted p-value).
#' Returns `mod_input` unchanged (with a warning) if `clusterProfiler` or
#' the annotation package is not installed.
#'
#' @param mod_input ModuleInput from [build_wgcna_modules()].
#' @param org_db Annotation package name. Default `"org.At.tair.db"`.
#' @param pval_cut FDR cutoff. Default `0.05`.
#' @return `mod_input` with `$module_meta$go_top` filled where enrichment
#'   was significant.
#' @export
annotate_go <- function(mod_input, org_db = "org.At.tair.db",
                        pval_cut = 0.05) {

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    warning("'clusterProfiler' is not installed. Skipping GO enrichment. ",
            "Install with: BiocManager::install('clusterProfiler')")
    return(mod_input)
  }
  if (!requireNamespace(org_db, quietly = TRUE)) {
    warning("'", org_db, "' is not installed. Skipping GO enrichment. ",
            "Install with: BiocManager::install('", org_db, "')")
    return(mod_input)
  }

  gene_mod    <- mod_input$gene_module
  all_genes   <- gene_mod$gene_id
  unique_mods <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))

  db <- get(org_db, envir = asNamespace(org_db))

  for (m in unique_mods) {
    mod_genes <- gene_mod$gene_id[gene_mod$top_module == m]
    tryCatch({
      ego <- clusterProfiler::enrichGO(
        gene          = mod_genes,
        universe      = all_genes,
        OrgDb         = db,
        keyType       = "TAIR",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cut,
        readable      = FALSE
      )
      if (!is.null(ego) && nrow(ego@result) > 0L) {
        best <- ego@result[which.min(ego@result$p.adjust), ]
        ri   <- which(mod_input$module_meta$module_id == m)
        mod_input$module_meta$go_top[ri] <- sprintf(
          "%s %s (p=%.3f)", best$ID, best$Description, best$p.adjust
        )
      }
    }, error = function(e) {
      message("GO enrichment failed for module ", m, ": ", conditionMessage(e))
    })
  }

  mod_input
}

# ---------------------------------------------------------------------------
# annotate_tfs
# ---------------------------------------------------------------------------

#' Intersect module genes with TF list
#'
#' Reads the lab TF metadata file (`Athaliana_motifs_metadata.tsv`, 673 TFs;
#' `motif_id` column = AT-ID) and fills `mod_input$module_tfs`. Use this
#' file — **NOT** PlantTFDB, which failed in the CZL run.
#'
#' @param mod_input ModuleInput from [build_wgcna_modules()].
#' @param tf_metadata_path Path to `Athaliana_motifs_metadata.tsv`.
#' @return `mod_input` with `$module_tfs` filled.
#' @export
annotate_tfs <- function(mod_input, tf_metadata_path) {
  if (!file.exists(tf_metadata_path))
    stop("TF metadata file not found: ", tf_metadata_path)

  tf_meta <- read.table(tf_metadata_path, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, quote = "")

  if (!("motif_id" %in% names(tf_meta)))
    stop("TF metadata must have a 'motif_id' column (AT-ID).")

  # Flexible column resolution
  sym_col <- if ("gene_symbol" %in% names(tf_meta)) "gene_symbol"
             else if ("symbol" %in% names(tf_meta)) "symbol"
             else NA_character_

  fam_col <- if ("tf_family" %in% names(tf_meta)) "tf_family"
             else if ("family" %in% names(tf_meta)) "family"
             else if ("class" %in% names(tf_meta)) "class"
             else NA_character_

  gene_mod    <- mod_input$gene_module
  unique_mods <- sort(unique(gene_mod$top_module[gene_mod$top_module > 0L]))

  tf_list <- lapply(unique_mods, function(m) {
    mod_genes <- gene_mod$gene_id[gene_mod$top_module == m]
    tf_in     <- tf_meta[tf_meta$motif_id %in% mod_genes, , drop = FALSE]
    if (nrow(tf_in) == 0L) return(NULL)
    data.frame(
      module_id   = as.integer(m),
      gene_id     = tf_in$motif_id,
      gene_symbol = if (!is.na(sym_col)) tf_in[[sym_col]] else NA_character_,
      tf_family   = if (!is.na(fam_col)) tf_in[[fam_col]] else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, Filter(Negate(is.null), tf_list))
  mod_input$module_tfs <- if (!is.null(result)) result else {
    data.frame(module_id = integer(), gene_id = character(),
               gene_symbol = character(), tf_family = character(),
               stringsAsFactors = FALSE)
  }

  mod_input
}
