## Profile GGM modules by condition pattern (FLAG-12)
## For each of the 4 official module sets, compute per-module condition-pattern
## distributions. Also produces BON3 (AT1G08860) condition-pattern table.
##
## Reads:  results/pathogen_multiome/robustness/pair_condition_patterns.csv
##         results/pathogen_multiome/official_modules/{set}/gene_module.csv
##
## Writes: results/pathogen_multiome/official_modules/{set}/module_condition_patterns.csv
##         results/pathogen_multiome/official_modules/all_modules_condition_patterns.csv
##         results/pathogen_multiome/robustness/BON3_condition_patterns.csv

suppressPackageStartupMessages(library(CoexprArabidopsis))

DATASET_ID  <- "pathogen_multiome"
RESULTS_DIR <- file.path("results", DATASET_ID)
ROB_DIR     <- file.path(RESULTS_DIR, "robustness")
MOD_DIR     <- file.path(RESULTS_DIR, "official_modules")
SETS        <- c("large_wgcna", "large_louvain", "small_wgcna", "small_louvain")

# Named pattern labels in display order (for stable column ordering)
NAMED_LABELS <- c("constitutive_all", "pan_pathogen", "ETI_shared",
                  "single_Mock", "single_DC3000", "single_AvrRpt2",
                  "single_AvrRpm1", "none")

# ---------------------------------------------------------------------------
message("\n=== Loading pair_condition_patterns.csv ===")
t0 <- proc.time()
cp <- read.csv(file.path(ROB_DIR, "pair_condition_patterns.csv"),
               stringsAsFactors = FALSE)
message(sprintf("Loaded %d pairs in %.1f s", nrow(cp), (proc.time()-t0)[["elapsed"]]))

# Pre-compute all unique pattern labels for column-ordering
all_labels  <- sort(unique(cp$pattern_label))
mixed_lbls  <- sort(grep("^mixed_", all_labels, value = TRUE))
label_order <- c(NAMED_LABELS[NAMED_LABELS %in% all_labels], mixed_lbls)

# ---------------------------------------------------------------------------
message("\n=== Processing module sets ===")

all_rows <- list()

for (set_name in SETS) {
  message("\n--- ", set_name, " ---")

  gm_path <- file.path(MOD_DIR, set_name, "gene_module.csv")
  gm <- read.csv(gm_path, stringsAsFactors = FALSE)

  # Exclude grey module (module 0 = unassigned genes)
  gm <- gm[gm$top_module != 0L, ]
  message(sprintf("  %d genes in %d non-grey modules", nrow(gm),
                  length(unique(gm$top_module))))

  gene_to_mod <- setNames(gm$top_module, gm$gene_id)

  # Add module assignments to all pairs; keep only intra-module pairs
  mod_A <- gene_to_mod[cp$gene_id_A]
  mod_B <- gene_to_mod[cp$gene_id_B]
  keep  <- !is.na(mod_A) & !is.na(mod_B) & (mod_A == mod_B)
  intra <- cp[keep, ]
  intra$module_id <- as.integer(mod_A[keep])

  message(sprintf("  %d intra-module pairs across %d modules",
                  nrow(intra), length(unique(intra$module_id))))

  if (nrow(intra) == 0L) {
    warning(set_name, ": no intra-module pairs found; skipping.")
    next
  }

  # Per-module statistics
  mods <- sort(unique(intra$module_id))
  rows <- lapply(mods, function(mid) {
    sub <- intra[intra$module_id == mid, ]
    n   <- nrow(sub)

    pat_tbl   <- table(sub$pattern_label)
    pat_frac  <- pat_tbl / n
    dom_pat   <- names(which.max(pat_frac))
    dom_frac  <- as.numeric(max(pat_frac))

    frac_cols <- setNames(
      vapply(label_order, function(lbl) {
        if (lbl %in% names(pat_frac)) as.numeric(pat_frac[[lbl]]) else 0.0
      }, numeric(1)),
      paste0("frac_", label_order)
    )

    r <- c(
      list(
        set           = set_name,
        module_id     = mid,
        n_intra_edges = n,
        dominant_pattern  = dom_pat,
        dominant_fraction = dom_frac,
        w_Mock    = mean(sub$w_Mock,    na.rm = TRUE),
        w_DC3000  = mean(sub$w_DC3000,  na.rm = TRUE),
        w_AvrRpt2 = mean(sub$w_AvrRpt2, na.rm = TRUE),
        w_AvrRpm1 = mean(sub$w_AvrRpm1, na.rm = TRUE),
        module_specificity_index   = mean(sub$specificity_index, na.rm = TRUE),
        n_conditions_active_mean   = mean(sub$n_conditions_active, na.rm = TRUE)
      ),
      as.list(frac_cols)
    )
    r
  })

  set_df <- do.call(rbind, lapply(rows, as.data.frame,
                                  stringsAsFactors = FALSE))

  out_path <- file.path(MOD_DIR, set_name, "module_condition_patterns.csv")
  write.csv(set_df, out_path, row.names = FALSE)
  message(sprintf("  Saved: %s", out_path))

  all_rows[[set_name]] <- set_df
}

# ---------------------------------------------------------------------------
message("\n=== Combined cross-set table ===")

all_df    <- do.call(rbind, all_rows)
combined_path <- file.path(MOD_DIR, "all_modules_condition_patterns.csv")
write.csv(all_df, combined_path, row.names = FALSE)
message("Saved: ", nrow(all_df), " module rows -> ", combined_path)

# Quick cross-set dominant-pattern summary
cat("\nDominant pattern counts per set:\n")
for (s in SETS) {
  sub <- all_df[all_df$set == s, ]
  tbl <- sort(table(sub$dominant_pattern), decreasing = TRUE)
  cat(sprintf("  %s (%d modules):\n", s, nrow(sub)))
  for (nm in names(tbl)) cat(sprintf("    %-24s %d\n", nm, tbl[[nm]]))
}

# ---------------------------------------------------------------------------
message("\n=== Part 3: BON3 (AT1G08860) condition patterns ===")

BON3_ID <- "AT1G08860"
bon3 <- cp[cp$gene_id_A == BON3_ID | cp$gene_id_B == BON3_ID, ]
message(sprintf("BON3 pairs found: %d", nrow(bon3)))

bon3_path <- file.path(ROB_DIR, "BON3_condition_patterns.csv")
write.csv(bon3, bon3_path, row.names = FALSE)
message("Saved: ", bon3_path)

cat("\nBON3 pattern summary:\n")
bon3_summary <- as.data.frame(table(pattern_label = bon3$pattern_label),
                              stringsAsFactors = FALSE)
bon3_summary <- bon3_summary[order(-bon3_summary$Freq), ]
print(bon3_summary)

message("\nDone.")
