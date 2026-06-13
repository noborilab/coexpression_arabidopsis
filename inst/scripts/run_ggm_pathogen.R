# inst/scripts/run_ggm_pathogen.R
# SingleCellGGM rerun — pathogen multiome, per condition
# Replaces the casual pooled run (FLAG-05, FLAG-06 in docs/PIPELINE_FLAGS.md)
#
# NOTE: this script constructs the InputBundle directly rather than calling
# load_seurat(), because the "condition" stratum must be derived from sample2
# (the column does not exist in the Seurat object and cannot be injected
# through the load_seurat() path argument interface without re-saving the .rds).
# All gene-ID mapping and min_cells logic mirrors what load_seurat() does.
#
# Usage — timing probe (Mock only, then report back):
#   Rscript inst/scripts/run_ggm_pathogen.R
#
# Usage — full run after timing confirmed:
#   Set CONDITIONS below to all 4, then:
#   nohup Rscript inst/scripts/run_ggm_pathogen.R > logs/ggm_rerun.log 2>&1 &

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(Matrix)
  library(jsonlite)
})

# --------------------------------------------------------------------------
# CONFIG — change CONDITIONS to c("Mock","DC3000","AvrRpt2","AvrRpm1") for
#           the full run after the timing probe confirms runtime is acceptable.
# --------------------------------------------------------------------------

CONDITIONS <- c("Mock", "DC3000", "AvrRpt2", "AvrRpm1")   # full run

SEURAT_PATH <- "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects/SA_PTI_ETI_single_cell/SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/motifFixed/combined_filtered.rds"

GFF3_PATH <- "/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/At_reference/Arabidopsis_thaliana.TAIR10.52.gff3"

OUTPUT_DIR <- "output_per_condition"
ASSAY      <- "RNA"
SLOT       <- "data"       # log-normalised layer
MIN_CELLS  <- 10L          # per-condition detection cutoff (matches coex_cutoff)
SEED       <- 98L
AT_ID_RE   <- "^AT[1-5MC]G[0-9]{5}$"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create("logs",     showWarnings = FALSE)

t0_total <- proc.time()[["elapsed"]]

# --------------------------------------------------------------------------
# 1. Load Seurat object
# --------------------------------------------------------------------------
message("[1/5] Loading Seurat object...")
obj <- readRDS(SEURAT_PATH)
message("  Loaded: ", nrow(obj), " genes x ", ncol(obj), " cells")

# --------------------------------------------------------------------------
# 2. Derive 'condition' column from sample2
#    sample2 values: "00_Mock" | "DC3000_04h" | "AvrRpt2_09h" | "AvrRpm1_24h"
#    → "Mock"        | "DC3000"               | "AvrRpt2"      | "AvrRpm1"
# --------------------------------------------------------------------------
message("[2/5] Deriving condition column from sample2...")
s2   <- as.character(obj@meta.data$sample2)
cond <- sub("_[0-9]{2}h.*$", "", sub("^00_", "", s2))
message("  Condition table:")
print(table(cond))

# --------------------------------------------------------------------------
# 3. Extract RNA data layer (SeuratObject v4/v5 agnostic)
# --------------------------------------------------------------------------
message("[3/5] Extracting RNA '", SLOT, "' layer...")
so_ver <- utils::packageVersion("SeuratObject")
counts <- if (so_ver >= "5.0.0") {
  SeuratObject::GetAssayData(obj, assay = ASSAY, layer = SLOT)
} else {
  SeuratObject::GetAssayData(obj, assay = ASSAY, slot  = SLOT)
}
message("  Matrix: ", nrow(counts), " x ", ncol(counts))

# --------------------------------------------------------------------------
# 4. Build symbol -> AT-ID map from TAIR10 GFF3 (base R, no rtracklayer)
#    Parse gene lines: ID=gene:ATID;Name=SYMBOL;...;gene_id=ATID;...
# --------------------------------------------------------------------------
message("[4/5] Building symbol -> AT-ID map from GFF3...")
gff_raw   <- readLines(GFF3_PATH)
gene_lines <- gff_raw[grepl("\tgene\t", gff_raw, fixed = TRUE)]
rm(gff_raw)
gc()
message("  Gene records in GFF3: ", length(gene_lines))

extract_attr <- function(lines, key) {
  pat <- paste0("(?:^|;)", key, "=([^;\\t]+)")
  m   <- regmatches(lines, regexpr(pat, lines, perl = TRUE))
  sub(paste0(".*", key, "="), "", m)
}

gff_gene_id <- extract_attr(gene_lines, "gene_id")
gff_name    <- extract_attr(gene_lines, "Name")

valid_at    <- grepl(AT_ID_RE, gff_gene_id)
symbol_map  <- data.frame(
  gene_symbol = gff_name[valid_at],
  gene_id     = gff_gene_id[valid_at],
  stringsAsFactors = FALSE
)
symbol_map <- symbol_map[!duplicated(symbol_map$gene_symbol), ]
message("  symbol_map: ", nrow(symbol_map), " unique symbol -> AT-ID pairs")

# Classify rownames
rn           <- rownames(counts)
is_at_id     <- grepl(AT_ID_RE, rn)
sym_rn       <- rn[!is_at_id]
map_idx      <- match(sym_rn, symbol_map$gene_symbol)
n_mapped     <- sum(!is.na(map_idx))
n_unmapped   <- sum(is.na(map_idx))

message("  Rowname gene ID breakdown:")
message("    Already AT-IDs:          ", sum(is_at_id))
message("    Symbols to map:          ", length(sym_rn))
message("    Mapped successfully:     ", n_mapped)
message("    Unmapped (dropped):      ", n_unmapped)
if (n_unmapped > 0 && n_unmapped <= 20) {
  message("    Unmapped: ", paste(sym_rn[is.na(map_idx)], collapse = ", "))
} else if (n_unmapped > 20) {
  message("    First 20 unmapped: ",
          paste(head(sym_rn[is.na(map_idx)], 20), collapse = ", "))
}

# Build final gene_id and gene_symbol vectors; drop unmapped.
# Use explicit integer indexing — double-subset assignment (x[a][b] <- v) does
# NOT write back to x in R, and OR-ing vectors of different lengths triggers
# silent recycling. Both are avoided here.
gene_id_all <- rn
sym_all     <- rep(NA_character_, length(rn))

sym_row_idx    <- which(!is_at_id)           # integer indices of symbol rows
mapped_sub_idx <- which(!is.na(map_idx))      # which of those were mapped

gene_id_all[sym_row_idx[mapped_sub_idx]] <-
  symbol_map$gene_id[map_idx[mapped_sub_idx]]
sym_all[sym_row_idx] <- sym_rn

keep_genes <- logical(length(rn))
keep_genes[is_at_id]                      <- TRUE   # all AT-IDs pass
keep_genes[sym_row_idx[mapped_sub_idx]]   <- TRUE   # mapped symbols pass

counts      <- counts[keep_genes, , drop = FALSE]
gene_id_all <- gene_id_all[keep_genes]
sym_all     <- sym_all[keep_genes]
rownames(counts) <- gene_id_all

gene_meta <- data.frame(
  gene_id     = gene_id_all,
  gene_symbol = sym_all,
  stringsAsFactors = FALSE
)

message("  Genes after mapping: ", nrow(counts))

# --------------------------------------------------------------------------
# 5. Cell metadata
# --------------------------------------------------------------------------
cell_meta          <- as.data.frame(obj@meta.data, stringsAsFactors = FALSE)
cell_meta$condition <- cond
cell_meta$cell_id  <- colnames(obj)
rm(obj); gc()

# --------------------------------------------------------------------------
# 6. Per-condition GGM
# --------------------------------------------------------------------------
message("[5/5] Running SingleCellGGM per condition...\n")

results_summary <- list()

for (this_cond in CONDITIONS) {
  message(strrep("=", 50))
  message("Condition: ", this_cond)
  message(strrep("=", 50))
  t0_cond <- proc.time()[["elapsed"]]

  keep_cells <- cell_meta$condition == this_cond
  if (sum(keep_cells) == 0) {
    warning("No cells for condition '", this_cond, "'; skipping.")
    next
  }

  counts_sub    <- counts[, keep_cells, drop = FALSE]
  meta_sub      <- cell_meta[keep_cells, , drop = FALSE]
  gene_meta_sub <- gene_meta  # same gene universe; GGM will filter per-stratum

  message("  Cells in stratum: ", sum(keep_cells))
  message("  Genes passed to GGM (pre-stratum filter): ", nrow(counts_sub))

  bundle <- list(
    counts       = counts_sub,
    cell_meta    = meta_sub,
    gene_meta    = gene_meta_sub,
    stratum_spec = list(variable = "condition", levels = this_cond),
    dataset_id   = "pathogen_multiome"
  )

  res <- estimate_singlecellggm(
    bundle        = bundle,
    n_iter        = NULL,
    subsample     = 2000L,
    pcor_cutoff   = 0.02,
    coex_cutoff   = MIN_CELLS,
    keep_negative = FALSE,
    seed          = SEED
  )

  nr       <- res[[this_cond]]
  t1_cond  <- proc.time()[["elapsed"]]
  wall_min <- round((t1_cond - t0_cond) / 60, 1)

  message("  n_iter:    ", nr$params$n_iter)
  message("  n_genes:   ", nr$params$n_genes, "  (post stratum filter)")
  message("  n_cells:   ", nr$params$n_cells)
  message("  n_edges:   ", nrow(nr$edge_table))
  message("  Wall time: ", wall_min, " min")

  out_dir <- file.path(OUTPUT_DIR, this_cond)
  dir.create(out_dir, showWarnings = FALSE)

  write.csv(nr$edge_table, file.path(out_dir, "edge_table.csv"), row.names = FALSE)

  params_out <- list(
    stratum_id              = this_cond,
    mode                    = "singlecellggm",
    n_cells                 = nr$params$n_cells,
    n_genes                 = nr$params$n_genes,
    n_unmapped_genes_dropped = n_unmapped,
    n_iter                  = nr$params$n_iter,
    subsample               = nr$params$subsample,
    aggregation             = "min_abs_pcor",
    pcor_cutoff             = nr$params$pcor_cutoff,
    coex_cutoff             = nr$params$coex_cutoff,
    keep_negative           = nr$params$keep_negative,
    seed                    = nr$params$seed,
    wall_time_min           = wall_min,
    timestamp               = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  jsonlite::write_json(params_out, file.path(out_dir, "params.json"),
                       pretty = TRUE, auto_unbox = TRUE)

  results_summary[[this_cond]] <- list(
    n_cells  = nr$params$n_cells,
    n_genes  = nr$params$n_genes,
    n_iter   = nr$params$n_iter,
    n_edges  = nrow(nr$edge_table),
    wall_min = wall_min
  )

  message("  Saved -> ", out_dir, "/edge_table.csv + params.json")
}

message("\n=== Run summary ===")
for (nm in names(results_summary)) {
  s <- results_summary[[nm]]
  message(sprintf("  %-10s  cells=%d  genes=%d  iter=%d  edges=%d  time=%.1f min",
                  nm, s$n_cells, s$n_genes, s$n_iter, s$n_edges, s$wall_min))
}
t1_total <- proc.time()[["elapsed"]]
message("Total elapsed: ", round((t1_total - t0_total) / 60, 1), " min")
message("Output dir: ", normalizePath(OUTPUT_DIR, mustWork = FALSE))
