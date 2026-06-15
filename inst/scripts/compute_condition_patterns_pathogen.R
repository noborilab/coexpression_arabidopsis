## Compute per-pair condition patterns for the pathogen multiome GGM
## Part 1 of FLAG-12: characterize_condition_pattern() on all 1.4M pairs.
## Output: results/pathogen_multiome/robustness/pair_condition_patterns.csv
##         results/pathogen_multiome/robustness/pattern_counts.csv

suppressPackageStartupMessages(library(CoexprArabidopsis))

DATASET_ID  <- "pathogen_multiome"
RESULTS_DIR <- file.path("results", DATASET_ID)
ROB_DIR     <- file.path(RESULTS_DIR, "robustness")
CONDITIONS  <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")

# Condition-pattern labels for this dataset. Supply to characterize_condition_pattern().
# For a new dataset, replace with labels appropriate to your condition set.
# Patterns not listed here fall back to the generic "pattern_<bits>" label.
PATTERN_LABELS <- c(
  "0000" = "none",           "1111" = "constitutive_all",
  "1000" = "single_Mock",    "0100" = "single_DC3000",
  "0010" = "single_AvrRpt2", "0001" = "single_AvrRpm1",
  "0111" = "pan_pathogen",   "0011" = "ETI_shared"
)

# ---------------------------------------------------------------------------
message("\n=== Loading inputs ===")

rob <- readRDS(file.path(ROB_DIR, "robustness_result.rds"))
message("rob: ", nrow(rob$pair_scores), " pairs")

network_list <- load_network_results(
  "output_per_condition",
  strata = CONDITIONS,
  mode   = "singlecellggm"
)
message("network_list: ", length(network_list), " conditions")

# ---------------------------------------------------------------------------
message("\n=== characterize_condition_pattern() ===")
t0 <- proc.time()

cp <- characterize_condition_pattern(rob, network_list, condition_order = CONDITIONS,
                                     pattern_labels = PATTERN_LABELS)

elapsed <- (proc.time() - t0)[["elapsed"]]
message(sprintf("Done in %.1f s — %d pairs x %d conditions", elapsed, nrow(cp), 4L))

# ---------------------------------------------------------------------------
message("\n=== Saving pair_condition_patterns.csv ===")

out_pairs <- file.path(ROB_DIR, "pair_condition_patterns.csv")
write.csv(cp, out_pairs, row.names = FALSE)
message("Saved: ", nrow(cp), " rows -> ", out_pairs)

# ---------------------------------------------------------------------------
message("\n=== Saving pattern_counts.csv ===")

pat_counts <- as.data.frame(table(
  pattern       = cp$pattern,
  pattern_label = cp$pattern_label
), stringsAsFactors = FALSE)
pat_counts <- pat_counts[pat_counts$Freq > 0L, ]
pat_counts <- pat_counts[order(-pat_counts$Freq), ]
names(pat_counts)[names(pat_counts) == "Freq"] <- "n_pairs"

out_counts <- file.path(ROB_DIR, "pattern_counts.csv")
write.csv(pat_counts, out_counts, row.names = FALSE)
message("Saved: ", nrow(pat_counts), " pattern entries -> ", out_counts)

# ---------------------------------------------------------------------------
message("\n=== Pattern summary ===")
print(pat_counts[, c("pattern_label", "pattern", "n_pairs")])

message("\nDone.")
