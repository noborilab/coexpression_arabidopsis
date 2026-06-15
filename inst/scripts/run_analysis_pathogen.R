#!/usr/bin/env Rscript
# run_analysis_pathogen.R
# End-to-end analysis: pathogen multiome (SingleCellGGM, per-condition)
# Steps: load → robustness → threshold sweep → WGCNA → annotate → save →
#        granularity sweep → cross-condition comparison → dev atlas probe → summary
#
# Run from project root:
#   nohup Rscript inst/scripts/run_analysis_pathogen.R > logs/analysis_pathogen.log 2>&1 &

suppressPackageStartupMessages({
  tryCatch(
    library(CoexprArabidopsis),
    error = function(e) {
      message("Package not installed; falling back to devtools::load_all()")
      devtools::load_all(".", quiet = TRUE)
    }
  )
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATASET_ID   <- "pathogen_multiome"
GGM_DIR      <- "output_per_condition"
RESULTS_DIR  <- file.path("results", DATASET_ID)
CONDITIONS   <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")
TF_META_PATH <- "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/TSL/from_Ben/for_tatsuya/data/motifs-2026/Athaliana_motifs_metadata.tsv"

for (d in c(file.path(RESULTS_DIR, "robustness"),
            file.path(RESULTS_DIR, "modules"),
            file.path(RESULTS_DIR, "granularity_sweep"),
            file.path(RESULTS_DIR, "condition_comparison"),
            "logs")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Tracking: completed steps for summary even if later steps fail
completed_steps <- character(0)

# ---------------------------------------------------------------------------
# Step 1 — load NetworkResults
# ---------------------------------------------------------------------------

message("\n=== Step 1: Load NetworkResults ===")
network_list <- tryCatch({
  nl <- load_network_results(GGM_DIR, strata = CONDITIONS, mode = "singlecellggm")
  for (nm in names(nl))
    message("  ", nm, ": ", nrow(nl[[nm]]$edge_table), " edges")
  nl
}, error = function(e) {
  message("FAILED step 1: ", conditionMessage(e))
  NULL
})

if (is.null(network_list)) stop("Cannot proceed without network_list.")

completed_steps <- c(completed_steps, "1_load")

# ---------------------------------------------------------------------------
# Step 2 — compute robustness
# ---------------------------------------------------------------------------

message("\n=== Step 2: Robustness ===")
rob <- tryCatch({
  r <- compute_robustness(network_list, k = 1.64, weight_cap = 30)
  save_robustness(r, file.path(RESULTS_DIR, "robustness"))
  message("Total pairs: ", nrow(r$pair_scores))
  print(quantile(r$pair_scores$R_score,
                 probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)))
  r
}, error = function(e) {
  message("FAILED step 2: ", conditionMessage(e))
  NULL
})

if (is.null(rob)) stop("Cannot proceed without robustness result.")

completed_steps <- c(completed_steps, "2_robustness")

# ---------------------------------------------------------------------------
# Step 3 — R_score threshold sweep
# ---------------------------------------------------------------------------

message("\n=== Step 3: R_score threshold sweep ===")
sweep_rscore <- NULL
R_SCORE_MIN  <- 0.5   # default

tryCatch({
  thresholds <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
  sweep_rscore <- data.frame(
    threshold = thresholds,
    n_edges   = sapply(thresholds, function(t) sum(rob$pair_scores$R_score >= t)),
    n_genes   = sapply(thresholds, function(t) {
      ps <- rob$pair_scores[rob$pair_scores$R_score >= t, ]
      length(unique(c(ps$gene_id_A, ps$gene_id_B)))
    })
  )
  print(sweep_rscore)
  write.csv(sweep_rscore,
            file.path(RESULTS_DIR, "robustness", "threshold_sweep.csv"),
            row.names = FALSE)

  ok <- sweep_rscore$n_genes >= 5000 & sweep_rscore$n_edges >= 50000
  R_SCORE_MIN <- if (any(ok)) min(sweep_rscore$threshold[ok]) else 0.5
  message("Auto-selected R_score threshold: ", R_SCORE_MIN)
  completed_steps <- c(completed_steps, "3_threshold_sweep")
}, error = function(e) {
  message("FAILED step 3: ", conditionMessage(e))
  message("Using default R_SCORE_MIN = ", R_SCORE_MIN)
})

# ---------------------------------------------------------------------------
# Step 4 — WGCNA primary modules
# ---------------------------------------------------------------------------

message("\n=== Step 4: WGCNA (primary) ===")
mod_input <- NULL
soft_power_used <- "auto"

tryCatch({
  # Capture messages to extract the auto-selected soft power
  msg_conn <- textConnection("msg_buffer", "w", local = TRUE)
  sink(msg_conn, type = "message")
  mi <- build_wgcna_modules(
    rob             = rob,
    network_list    = network_list,
    r_score_min     = R_SCORE_MIN,
    soft_power      = NULL,
    merge_cut       = 0.25,
    min_module_size = 30,
    sub_merge_cut   = 0.10
  )
  sink(NULL, type = "message")
  close(msg_conn)
  # Parse soft power from captured messages
  sp_line <- grep("Auto-selected soft_power", msg_buffer, value = TRUE)
  if (length(sp_line) > 0) {
    sp_val <- regmatches(sp_line[1], regexpr("[0-9]+", sp_line[1]))
    if (length(sp_val) > 0) soft_power_used <- sp_val
  }
  # Replay captured messages to the log
  for (m in msg_buffer) message(m)

  mod_input <- mi
  message("Top modules: ",
    length(unique(mod_input$gene_module$top_module[
      mod_input$gene_module$top_module != 0])))
  message("Genes assigned: ", sum(mod_input$gene_module$top_module != 0))
  completed_steps <- c(completed_steps, "4_wgcna_primary")
}, error = function(e) {
  # Ensure sink is always closed
  try(sink(NULL, type = "message"), silent = TRUE)
  message("FAILED step 4: ", conditionMessage(e))
})

if (is.null(mod_input)) stop("Cannot proceed without module assignments.")

# ---------------------------------------------------------------------------
# Step 5 — annotate primary modules
# ---------------------------------------------------------------------------

message("\n=== Step 5: Annotation ===")
tryCatch({
  mod_input <- annotate_context(mod_input, network_list, ref_condition = "Mock")
  mod_input <- annotate_go(mod_input, org_db = "org.At.tair.db", pval_cut = 0.05)

  if (file.exists(TF_META_PATH)) {
    mod_input <- annotate_tfs(mod_input, TF_META_PATH)
    message("TF entries: ", nrow(mod_input$module_tfs))
  } else {
    message("WARNING: TF file not found: ", TF_META_PATH)
  }

  if (length(network_list) > 1) {
    pres <- compute_preservation_fallback(
      mod_input     = mod_input,
      network_list2 = network_list[names(network_list) != "Mock"]
    )
    mod_input$module_meta <- merge(mod_input$module_meta, pres,
                                    by = "module_id", all.x = TRUE)
  }
  completed_steps <- c(completed_steps, "5_annotation")
}, error = function(e) {
  message("FAILED step 5: ", conditionMessage(e))
})

# ---------------------------------------------------------------------------
# Step 6 — save primary ModuleInput
# ---------------------------------------------------------------------------

message("\n=== Step 6: Save ===")
tryCatch({
  mdir <- file.path(RESULTS_DIR, "modules")
  write.csv(mod_input$gene_module,
            file.path(mdir, "gene_module.csv"),  row.names = FALSE)
  write.csv(mod_input$module_meta,
            file.path(mdir, "module_meta.csv"),  row.names = FALSE)
  write.csv(mod_input$module_hier,
            file.path(mdir, "module_hier.csv"),  row.names = FALSE)
  write.csv(mod_input$hub_genes,
            file.path(mdir, "hub_genes.csv"),    row.names = FALSE)
  write.csv(mod_input$module_tfs,
            file.path(mdir, "module_tfs.csv"),   row.names = FALSE)
  write.csv(as.data.frame(mod_input$eigengenes),
            file.path(mdir, "eigengenes.csv"),   row.names = TRUE)
  saveRDS(mod_input, file.path(mdir, "module_input.rds"))
  message("ModuleInput saved.")
  completed_steps <- c(completed_steps, "6_save_primary")
}, error = function(e) {
  message("FAILED step 6: ", conditionMessage(e))
})

# ---------------------------------------------------------------------------
# Step 7 — granularity sweep
# ---------------------------------------------------------------------------

message("\n=== Step 7: Granularity sweep ===")
sweep_df <- NULL

tryCatch({
  merge_cuts <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35)
  sweep_results <- list()

  t_sweep_start <- proc.time()[["elapsed"]]

  for (mc in merge_cuts) {
    for (sp in list(NULL, 4L, 6L, 8L)) {
      # Check elapsed time; trim grid if we're over 150 min total
      elapsed_min <- (proc.time()[["elapsed"]] - t_sweep_start) / 60
      if (elapsed_min > 150) {
        message("  Sweep running > 150 min — stopping grid early.")
        break
      }

      sp_label <- if (is.null(sp)) "auto" else as.character(sp)
      tag <- paste0("merge", mc, "_power", sp_label)
      message("  Running: ", tag)
      tryCatch({
        mi_sweep <- build_wgcna_modules(
          rob             = rob,
          network_list    = network_list,
          r_score_min     = R_SCORE_MIN,
          soft_power      = sp,
          merge_cut       = mc,
          min_module_size = 30,
          sub_merge_cut   = NULL
        )
        n_top <- length(unique(
          mi_sweep$gene_module$top_module[
            mi_sweep$gene_module$top_module != 0]))
        n_assigned <- sum(mi_sweep$gene_module$top_module != 0)
        n_grey     <- sum(mi_sweep$gene_module$top_module == 0)
        sweep_results[[tag]] <- data.frame(
          merge_cut  = mc,
          soft_power = sp_label,
          n_modules  = n_top,
          n_assigned = n_assigned,
          n_grey     = n_grey,
          pct_grey   = round(100 * n_grey / (n_assigned + n_grey), 1)
        )
      }, error = function(e) {
        message("  FAILED: ", tag, " — ", conditionMessage(e))
        sweep_results[[tag]] <<- data.frame(
          merge_cut = mc, soft_power = sp_label,
          n_modules = NA, n_assigned = NA, n_grey = NA, pct_grey = NA)
      })
    }
    # Check elapsed again at outer loop boundary
    elapsed_min <- (proc.time()[["elapsed"]] - t_sweep_start) / 60
    if (elapsed_min > 150) break
  }

  sweep_df <- do.call(rbind, sweep_results)
  rownames(sweep_df) <- NULL
  print(sweep_df)
  write.csv(sweep_df,
            file.path(RESULTS_DIR, "granularity_sweep", "sweep_results.csv"),
            row.names = FALSE)
  message("Granularity sweep saved.")
  completed_steps <- c(completed_steps, "7_granularity_sweep")
}, error = function(e) {
  message("FAILED step 7: ", conditionMessage(e))
})

# ---------------------------------------------------------------------------
# Step 8 — cross-condition module comparison
# ---------------------------------------------------------------------------

message("\n=== Step 8: Cross-condition comparison ===")
cond_df <- NULL

tryCatch({
  gene_mod <- mod_input$gene_module
  modules  <- unique(gene_mod$top_module[gene_mod$top_module != 0])

  cond_profiles <- lapply(modules, function(m) {
    genes_m <- gene_mod$gene_id[gene_mod$top_module == m]
    per_cond <- sapply(names(network_list), function(cond) {
      et    <- network_list[[cond]]$edge_table
      intra <- et[et$gene_id_A %in% genes_m & et$gene_id_B %in% genes_m, ]
      if (nrow(intra) == 0) return(NA_real_)
      mean(abs(intra$weight))
    })
    data.frame(module = m, t(per_cond), check.names = FALSE)
  })
  cond_df <- do.call(rbind, cond_profiles)

  cond_cols <- names(network_list)
  cond_df$top_condition <- apply(cond_df[, cond_cols], 1, function(x) {
    if (all(is.na(x))) return(NA_character_)
    cond_cols[which.max(x)]
  })

  for (cond in setdiff(cond_cols, "Mock")) {
    cond_df[[paste0("delta_", cond)]] <-
      cond_df[[cond]] - cond_df[["Mock"]]
  }

  print(cond_df)
  write.csv(cond_df,
            file.path(RESULTS_DIR, "condition_comparison",
                      "module_condition_profiles.csv"),
            row.names = FALSE)

  delta_cols <- paste0("delta_", setdiff(cond_cols, "Mock"))
  gained <- cond_df[apply(
    cond_df[, delta_cols, drop = FALSE],
    1, function(x) any(!is.na(x) & x > 0.005)), ]
  message("Modules gained in at least one pathogen condition: ", nrow(gained))
  write.csv(gained,
            file.path(RESULTS_DIR, "condition_comparison",
                      "modules_gained_in_pathogen.csv"),
            row.names = FALSE)
  completed_steps <- c(completed_steps, "8_condition_comparison")
}, error = function(e) {
  message("FAILED step 8: ", conditionMessage(e))
})

# ---------------------------------------------------------------------------
# Step 9 — dev atlas probe
# ---------------------------------------------------------------------------

message("\n=== Step 9: Dev atlas probe ===")
tryCatch({
  DEV_ATLAS_CANDIDATES <- c(
    "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects/sl_atlas",
    "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects/SL_atlas",
    "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/TSL/projects/sl_atlas",
    "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/TSL/projects/SL_atlas"
  )

  for (candidate in DEV_ATLAS_CANDIDATES) {
    if (dir.exists(candidate)) {
      message("Found dir: ", candidate)
      rds_files <- list.files(candidate, pattern = "\\.rds$",
                               recursive = TRUE, full.names = TRUE)
      if (length(rds_files) > 0) {
        sizes <- file.size(rds_files)
        message("  RDS files:")
        for (i in seq_along(rds_files)) {
          message("    ", rds_files[i], " (", round(sizes[i]/1e9, 2), " GB)")
        }
      } else {
        message("  No .rds files found")
      }
    }
  }

  message("\nSearching for pseudobulk count tables (dev atlas):")
  search_roots <- c(
    "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects",
    "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/TSL/projects"
  )
  for (root in search_roots) {
    if (!dir.exists(root)) next
    csvs <- list.files(root, pattern = "pseudobulk.*\\.csv$|.*pseudobulk\\.csv$",
                       recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(csvs) > 0) {
      message("  Found in ", root, ":")
      for (f in csvs) message("    ", f, " (", round(file.size(f)/1e6, 1), " MB)")
    }
  }
  message("Dev atlas probe complete.")
  completed_steps <- c(completed_steps, "9_dev_atlas_probe")
}, error = function(e) {
  message("FAILED step 9: ", conditionMessage(e))
})

# ---------------------------------------------------------------------------
# Step 10 — write ANALYSIS_SUMMARY.md
# ---------------------------------------------------------------------------

message("\n=== Step 10: Summary ===")
tryCatch({
  summary_path <- file.path(RESULTS_DIR, "ANALYSIS_SUMMARY.md")

  # Module meta columns actually present (they vary with annotation outcomes)
  mm_cols <- intersect(
    c("module_id","n_genes","top_organ_or_condition","go_top","zsummary","preservation_method"),
    colnames(mod_input$module_meta)
  )

  lines <- c(
    "# Pathogen Multiome — Analysis Summary",
    paste0("Date: ", Sys.time()),
    "",
    "## Run parameters",
    paste0("- Dataset: ", DATASET_ID),
    paste0("- R_score threshold: ", R_SCORE_MIN),
    paste0("- Soft power: ", soft_power_used),
    paste0("- merge_cut (top): 0.25 | sub_merge_cut: 0.10"),
    paste0("- min_module_size: 30"),
    "",
    "## Input",
    paste(sapply(names(network_list), function(nm)
      paste0("- ", nm, ": ", nrow(network_list[[nm]]$edge_table), " edges")),
      collapse = "\n"),
    "",
    "## Robustness",
    paste0("- Total pairs tested: ", nrow(rob$pair_scores)),
    paste0("- R_score quantiles:"),
    paste(capture.output(
      print(quantile(rob$pair_scores$R_score,
                     probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)))),
      collapse = "\n"),
    "",
    "## Threshold sweep",
    if (!is.null(sweep_rscore)) paste(capture.output(print(sweep_rscore)), collapse = "\n")
    else "(not computed)",
    "",
    "## Primary modules",
    paste(capture.output(print(
      mod_input$module_meta[, mm_cols, drop = FALSE]
    )), collapse = "\n"),
    "",
    "## Granularity sweep",
    if (!is.null(sweep_df)) paste(capture.output(print(sweep_df)), collapse = "\n")
    else "(not computed)",
    "",
    "## Cross-condition profiles",
    if (!is.null(cond_df)) paste(capture.output(print(cond_df)), collapse = "\n")
    else "(not computed)",
    "",
    "## Top 10 edges by R_score",
    paste(capture.output(print(
      head(rob$pair_scores[order(-rob$pair_scores$R_score),
           c("gene_id_A","gene_id_B","R_score","z_bar")], 10)
    )), collapse = "\n"),
    "",
    "## Hub genes (top 5 per module by kME)",
    paste(capture.output(print(
      do.call(rbind, lapply(split(mod_input$hub_genes,
                                  mod_input$hub_genes$module_id),
                            function(x) head(x[order(-x$kME),], 5)))
    )), collapse = "\n"),
    "",
    "## Completed steps",
    paste(completed_steps, collapse = ", ")
  )

  writeLines(lines, summary_path)
  message("Summary written to ", summary_path)
  completed_steps <- c(completed_steps, "10_summary")
}, error = function(e) {
  message("FAILED step 10: ", conditionMessage(e))
})

message("\n=== PIPELINE COMPLETE ===")
message("Completed steps: ", paste(completed_steps, collapse = ", "))
