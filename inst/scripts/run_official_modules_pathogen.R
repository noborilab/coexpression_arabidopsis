## Official module sets for pathogen_multiome
## Four sets: (large|small graph) x (WGCNA p1|Louvain)
## Run under nohup; no wall-time limit.

suppressPackageStartupMessages(library(CoexprArabidopsis))
suppressPackageStartupMessages(library(igraph))

t_global_start <- proc.time()

DATASET_ID  <- "pathogen_multiome"
RESULTS_DIR <- file.path("results", DATASET_ID)
OUT_DIR     <- file.path(RESULTS_DIR, "official_modules")
CONDITIONS  <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")

TF_META_PATH <- paste0(
  "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Nobori Lab (TSL) Team Folder/",
  "shared/datasets/from_Ben/for_tatsuya/data/motifs-2026/",
  "Athaliana_motifs_metadata.tsv"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load shared inputs once
# ---------------------------------------------------------------------------

message("\n=== Loading shared inputs ===")

rob <- readRDS(file.path(RESULTS_DIR, "robustness", "robustness_result.rds"))
message("rob: ", nrow(rob$pair_scores), " pairs")

network_list <- load_network_results(
  "output_per_condition",
  strata = CONDITIONS,
  mode   = "singlecellggm"
)
message("network_list: ", length(network_list), " conditions")

symbol_map <- read.csv(
  file.path(RESULTS_DIR, "symbol_map.csv"),
  stringsAsFactors = FALSE
)
message("symbol_map: ", nrow(symbol_map), " entries")

# Build named lookup for fast symbol resolution
sym_lookup <- setNames(symbol_map$gene_symbol, symbol_map$gene_id)

# ---------------------------------------------------------------------------
# Helper: join gene_symbol onto a data.frame that has a gene_id column
# ---------------------------------------------------------------------------

.add_symbol <- function(df, col = "gene_id") {
  if (!(col %in% names(df))) return(df)
  df$gene_symbol <- sym_lookup[df[[col]]]
  df
}

# ---------------------------------------------------------------------------
# Helper: build Louvain ModuleInput from filtered pair_scores
# ---------------------------------------------------------------------------

.build_louvain_modules <- function(rob, r_score_min) {

  ps <- rob$pair_scores
  ps <- ps[!is.na(ps$R_score) & ps$R_score >= r_score_min, , drop = FALSE]
  if (nrow(ps) == 0L)
    stop("No pairs remain after r_score_min filter = ", r_score_min)

  edges <- data.frame(
    from   = ps$gene_id_A,
    to     = ps$gene_id_B,
    weight = abs(tanh(ps$z_bar)),
    stringsAsFactors = FALSE
  )
  edges$weight <- pmin(edges$weight, 1.0)

  g  <- igraph::graph_from_data_frame(edges, directed = FALSE)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)

  memb    <- igraph::membership(cl)
  gene_ids <- names(memb)

  # Map small communities (< 30 genes) to module 0
  comm_sizes <- table(memb)
  small_comm <- as.integer(names(comm_sizes[comm_sizes < 30L]))
  top_labels <- as.integer(memb)
  top_labels[top_labels %in% small_comm] <- 0L

  # Re-label remaining modules as consecutive integers
  live_comms <- sort(unique(top_labels[top_labels > 0L]))
  relabel    <- setNames(seq_along(live_comms), as.character(live_comms))
  top_labels_final <- ifelse(top_labels == 0L, 0L,
                             as.integer(relabel[as.character(top_labels)]))

  # kME for Louvain: correlate each gene's adjacency row with
  # mean adjacency profile of its module members
  all_genes   <- sort(union(edges$from, edges$to))
  n_genes     <- length(all_genes)
  gene_idx    <- setNames(seq_along(all_genes), all_genes)

  A <- matrix(0.0, nrow = n_genes, ncol = n_genes,
              dimnames = list(all_genes, all_genes))
  i_A <- gene_idx[edges$from]
  i_B <- gene_idx[edges$to]
  A[cbind(i_A, i_B)] <- edges$weight
  A[cbind(i_B, i_A)] <- edges$weight

  # Reorder top_labels_final to match all_genes order
  tl_named <- setNames(top_labels_final, gene_ids)
  tl_ord   <- tl_named[all_genes]

  unique_mods <- sort(unique(tl_ord[tl_ord > 0L]))

  kme_vec <- rep(NA_real_, n_genes)
  for (m in unique_mods) {
    mod_idx  <- which(tl_ord == m)
    if (length(mod_idx) < 2L) next
    mod_mean <- colMeans(A[mod_idx, , drop = FALSE])
    for (i in mod_idx) {
      kme_vec[i] <- cor(A[i, ], mod_mean)
    }
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

  # Louvain has no hierarchy; return empty module_hier
  module_hier <- data.frame(
    sub_module = integer(),
    top_module = integer(),
    stringsAsFactors = FALSE
  )

  # hub_genes: top 20 per module by kME
  hub_list <- lapply(as.integer(names(mod_counts)), function(m) {
    sub_gm <- gene_module[!is.na(gene_module$kME) &
                          gene_module$top_module == m, , drop = FALSE]
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

  module_tfs <- data.frame(
    module_id = integer(), gene_id = character(),
    gene_symbol = character(), tf_family = character(),
    stringsAsFactors = FALSE
  )

  # eigengenes: mean adjacency profile per module (samples = modules themselves)
  # Use mean per-module expression proxy from adjacency
  ME_mat <- matrix(NA_real_, nrow = length(unique_mods), ncol = n_genes,
                   dimnames = list(paste0("ME", unique_mods), all_genes))
  for (m in unique_mods) {
    mod_idx <- which(tl_ord == m)
    ME_mat[paste0("ME", m), ] <- colMeans(A[mod_idx, , drop = FALSE])
  }
  # eigengenes slot: samples x modules (transpose)
  MEs <- as.data.frame(t(ME_mat))

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
# Helper: annotate and save one module set
# ---------------------------------------------------------------------------

.run_set <- function(set_name, mod_input, network_list, rob, tf_meta_path,
                     out_dir) {

  t0 <- proc.time()
  set_dir <- file.path(out_dir, set_name)
  dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n--- Annotating: ", set_name, " ---")

  # Add gene_symbol to gene_module and hub_genes
  mod_input$gene_module <- .add_symbol(mod_input$gene_module)
  mod_input$hub_genes   <- .add_symbol(mod_input$hub_genes)

  # 3a. Condition context
  mod_input <- tryCatch(
    annotate_context(mod_input, network_list, ref_condition = "Mock"),
    error = function(e) {
      message("annotate_context failed: ", conditionMessage(e)); mod_input
    }
  )

  # 3b. GO BP enrichment
  mod_input <- tryCatch(
    annotate_go(mod_input, org_db = "org.At.tair.db", pval_cut = 0.05),
    error = function(e) {
      message("annotate_go failed: ", conditionMessage(e)); mod_input
    }
  )

  # 3c. TF intersection
  if (file.exists(tf_meta_path)) {
    mod_input <- tryCatch(
      annotate_tfs(mod_input, tf_meta_path),
      error = function(e) {
        message("annotate_tfs failed: ", conditionMessage(e)); mod_input
      }
    )
  } else {
    message("TF file not found: ", tf_meta_path)
  }

  # 3d. Module preservation fallback (Mock as reference; test in pathogen conditions)
  test_nets <- network_list[names(network_list) != "Mock"]
  pres <- tryCatch(
    compute_preservation_fallback(mod_input, test_nets),
    error = function(e) {
      message("compute_preservation_fallback failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(pres)) {
    mod_input$module_meta <- merge(
      mod_input$module_meta, pres,
      by = "module_id", all.x = TRUE
    )
    # merge may create .x / .y duplicates; resolve
    if ("zsummary.x" %in% names(mod_input$module_meta)) {
      mod_input$module_meta$zsummary <- mod_input$module_meta$zsummary.y
      mod_input$module_meta$preservation_method <-
        mod_input$module_meta$preservation_method.y
      mod_input$module_meta$zsummary.x               <- NULL
      mod_input$module_meta$zsummary.y               <- NULL
      mod_input$module_meta$preservation_method.x    <- NULL
      mod_input$module_meta$preservation_method.y    <- NULL
    }
  }

  # 4. Save
  write.csv(mod_input$gene_module,
            file.path(set_dir, "gene_module.csv"),  row.names = FALSE)
  write.csv(mod_input$module_meta,
            file.path(set_dir, "module_meta.csv"),  row.names = FALSE)
  write.csv(mod_input$module_hier,
            file.path(set_dir, "module_hier.csv"),  row.names = FALSE)
  write.csv(mod_input$hub_genes,
            file.path(set_dir, "hub_genes.csv"),    row.names = FALSE)
  write.csv(mod_input$module_tfs,
            file.path(set_dir, "module_tfs.csv"),   row.names = FALSE)
  write.csv(as.data.frame(mod_input$eigengenes),
            file.path(set_dir, "eigengenes.csv"),   row.names = TRUE)
  saveRDS(mod_input, file.path(set_dir, "module_input.rds"))

  elapsed <- (proc.time() - t0)[["elapsed"]]
  message("Saved: ", set_name, " in ", round(elapsed / 60, 1), " min")

  mod_input
}

# ---------------------------------------------------------------------------
# Build and annotate all four official sets
# ---------------------------------------------------------------------------

results_list   <- list()
timing_list    <- list()
set_configs <- list(
  large_wgcna   = list(graph = "large", method = "wgcna",   r_min = 0.5),
  large_louvain  = list(graph = "large", method = "louvain", r_min = 0.5),
  small_wgcna   = list(graph = "small", method = "wgcna",   r_min = 0.6),
  small_louvain  = list(graph = "small", method = "louvain", r_min = 0.6)
)

for (set_name in names(set_configs)) {

  cfg <- set_configs[[set_name]]
  message("\n==============================")
  message("=== SET: ", set_name, " (", cfg$method, ", R>=", cfg$r_min, ") ===")
  message("==============================")
  t_set <- proc.time()

  tryCatch({

    # 1. Build modules
    mod_input <- if (cfg$method == "wgcna") {
      build_wgcna_modules(
        rob             = rob,
        network_list    = network_list,
        r_score_min     = cfg$r_min,
        soft_power      = 1L,
        merge_cut       = 0.25,
        min_module_size = 30,
        sub_merge_cut   = 0.10
      )
    } else {
      .build_louvain_modules(rob, r_score_min = cfg$r_min)
    }

    # Add method/graph metadata fields to mod_input
    mod_input$method          <- if (cfg$method == "wgcna") "wgcna_p1" else "louvain"
    mod_input$graph           <- cfg$graph
    mod_input$r_score_threshold <- cfg$r_min

    # Quick stats before annotation
    gm <- mod_input$gene_module
    n_modules  <- length(unique(gm$top_module[gm$top_module > 0L]))
    n_assigned <- sum(gm$top_module > 0L)
    n_grey     <- sum(is.na(gm$top_module) | gm$top_module == 0L)
    n_total    <- nrow(gm)
    pct_grey   <- round(100 * n_grey / n_total, 1)
    message(sprintf("  Modules: %d | Assigned: %d | Grey/unassigned: %d (%.1f%%)",
                    n_modules, n_assigned, n_grey, pct_grey))

    # 2–4. Annotate + save
    mod_input <- .run_set(
      set_name     = set_name,
      mod_input    = mod_input,
      network_list = network_list,
      rob          = rob,
      tf_meta_path = TF_META_PATH,
      out_dir      = OUT_DIR
    )

    elapsed_set <- (proc.time() - t_set)[["elapsed"]]
    timing_list[[set_name]] <- elapsed_set
    results_list[[set_name]] <- mod_input

    message(sprintf("  TF entries: %d", nrow(mod_input$module_tfs)))

  }, error = function(e) {
    message("ERROR in set '", set_name, "': ", conditionMessage(e))
    elapsed_set <- (proc.time() - t_set)[["elapsed"]]
    timing_list[[set_name]] <<- elapsed_set
  })
}

# ---------------------------------------------------------------------------
# Cross-set gene assignment comparison
# ---------------------------------------------------------------------------

message("\n=== Cross-set gene assignment comparison ===")

if (length(results_list) == 4L) {
  tryCatch({
    gm_list <- lapply(names(results_list), function(sn) {
      gm <- results_list[[sn]]$gene_module[, c("gene_id", "top_module"), drop = FALSE]
      names(gm)[2] <- paste0("module_", sn)
      gm
    })

    cross <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE), gm_list)
    cross$gene_symbol <- sym_lookup[cross$gene_id]

    # Move gene_symbol to second column
    cross <- cross[, c("gene_id", "gene_symbol",
                       "module_large_wgcna", "module_large_louvain",
                       "module_small_wgcna", "module_small_louvain")]

    write.csv(cross,
              file.path(OUT_DIR, "cross_set_assignments.csv"),
              row.names = FALSE)

    # Summary: how many genes have the same top_module > 0 assignment in all 4?
    # (Louvain and WGCNA module IDs are not comparable, so "agreement" means
    # all four assign the gene to a non-grey module — a coverage metric.)
    all_assigned <- rowSums(!is.na(cross[, 3:6]) &
                            cross[, 3:6] > 0, na.rm = TRUE)
    message(sprintf(
      "cross_set_assignments.csv: %d genes | assigned in all 4 sets: %d | grey in at least one: %d",
      nrow(cross),
      sum(all_assigned == 4L, na.rm = TRUE),
      sum(all_assigned < 4L, na.rm = TRUE)
    ))
  }, error = function(e) {
    message("Cross-set comparison failed: ", conditionMessage(e))
  })
} else {
  message("Not all 4 sets completed; skipping cross-set comparison.")
}

# ---------------------------------------------------------------------------
# Completion report
# ---------------------------------------------------------------------------

t_total <- (proc.time() - t_global_start)[["elapsed"]]

message("\n========== COMPLETION REPORT ==========")
for (set_name in names(set_configs)) {
  mi <- results_list[[set_name]]
  if (is.null(mi)) {
    message(sprintf("%-18s: FAILED", set_name))
    next
  }
  gm         <- mi$gene_module
  n_modules  <- length(unique(gm$top_module[gm$top_module > 0L]))
  n_assigned <- sum(!is.na(gm$top_module) & gm$top_module > 0L)
  n_grey     <- sum(is.na(gm$top_module) | gm$top_module == 0L)
  n_total    <- nrow(gm)
  pct_grey   <- round(100 * n_grey / n_total, 1)
  n_tf       <- nrow(mi$module_tfs)
  elapsed    <- timing_list[[set_name]]
  top_go_sample <- if ("go_top" %in% names(mi$module_meta)) {
    go_vals <- mi$module_meta$go_top[!is.na(mi$module_meta$go_top)]
    if (length(go_vals) > 0) go_vals[1] else "none"
  } else "none"

  message(sprintf(
    "%-18s | modules=%d | assigned=%d | grey=%d (%.1f%%) | TF_entries=%d | wall=%.1fmin | top_GO_m1='%s'",
    set_name, n_modules, n_assigned, n_grey, pct_grey, n_tf,
    elapsed / 60, top_go_sample
  ))
}

tf_found <- file.exists(TF_META_PATH)
message(sprintf("\nTF file found: %s", tf_found))
message(sprintf("Total wall time: %.1f min", t_total / 60))
message("=== DONE ===")
