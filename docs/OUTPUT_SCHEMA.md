# Output Schema Reference

This document describes every data object the pipeline produces: what it
contains, where it comes from, where it goes, and how it is saved. A wet-lab
collaborator who knows basic R should be able to read this and understand what
file to open and what each column means.

---

## InputBundle

**One-line description:** Standardised input package that all core functions
receive. Seurat-specific logic ends here.

**Produced by:** `R/adapter_seurat.R`  
**Consumed by:** `R/estimate_pseudobulk.R`, `R/estimate_singlecellggm.R`  
**Saved to disk:** Not saved (transient in-memory object). The Seurat object it
is derived from is the authoritative source.

### Slots

| Slot | Type | Description |
|---|---|---|
| `$counts` | matrix (genes × cells) | Log-normalised RNA counts. Rows are genes (AT-IDs as rownames), columns are cells (barcodes as colnames). Genes passing the expression filter (detected in ≥ min_cells cells within the target subset) only. |
| `$cell_meta` | data.frame | One row per cell. Always contains `cell_id`, `stratum_var` (value of the condition/organ variable, e.g. "Mock"), and `group_var` (pseudobulk grouping label). Additional Seurat metadata columns passed through. |
| `$gene_meta` | data.frame | One row per gene in `$counts`. Columns: `gene_id` (AT-ID, Araport11), `gene_symbol` (gene symbol; NA if unmapped). |
| `$stratum_spec` | named list | `$variable` (character): column name in `cell_meta` used to stratify (e.g. "condition"). `$levels` (character vector): ordered condition/organ levels to process. |
| `$dataset_id` | character | Short identifier for the dataset, e.g. "pathogen_multiome". Propagated into all downstream objects. |

---

## NetworkResult

**One-line description:** A single co-expression network for one stratum level
(one condition or one organ), produced by either estimation mode.

**Produced by:** `R/estimate_pseudobulk.R` or `R/estimate_singlecellggm.R`
(each function returns a **list** of NetworkResult, one per stratum level)  
**Consumed by:** `R/robustness.R` (when enabled), `R/interpret.R` (directly
when robustness is skipped)  
**Saved to disk:**

- Edge table → `output_per_condition/{stratum_id}/edge_table.csv` (CSV)
- Parameters → `output_per_condition/{stratum_id}/params.json` (JSON sidecar)

### Slots

| Slot | Type | Description |
|---|---|---|
| `$edge_table` | data.frame | One row per gene pair with a non-zero edge. Columns: `gene_id_A` (AT-ID), `gene_id_B` (AT-ID), `weight` (pcor for GGM; Spearman r for pseudobulk). No gene symbols. |
| `$gene_ids` | character vector | AT-IDs of all genes in the network (i.e. all genes that passed the expression filter for this stratum). |
| `$stratum_id` | character | Which stratum this network represents, e.g. "Mock", "DC3000". |
| `$mode` | character | Either "pseudobulk" or "singlecellggm". |
| `$params` | named list | All parameters used to produce this network: `pcor_cutoff`, `n_iter`, `subsample`, `aggregation`, `min_cells`, `seed` for GGM; `min_pseudobulk_samples` for pseudobulk. |
| `$timestamp` | POSIXct | When this network was produced. |

### Notes on `$edge_table`

- For SingleCellGGM: `weight` is the minimum absolute partial correlation
  across iterations (the conservative, defining aggregation of the method).
  Only pairs with `weight ≥ pcor_cutoff` AND detected together in ≥ `min_cells`
  cells are included.
- For pseudobulk: `weight` is Spearman r (computed via rank-transform-then-Pearson
  on the pseudobulk matrix).

---

## RobustnessResult

**One-line description:** Cross-stratum reproducibility statistics for every
gene pair tested, across all stratum levels of one dataset.

**Produced by:** `R/robustness.R`  
**Consumed by:** `R/interpret.R` (to assemble ModuleInput)  
**Saved to disk:**

- `results/{dataset_id}/robustness/pair_scores_full.csv` — **all** tested pairs (CSV)
- `results/{dataset_id}/robustness/robustness_result.rds` — full R object (RDS)

**IMPORTANT:** The CSV must contain all tested pairs, not just high-scoring
ones. In the prior CZL run, only R ≥ 0.7 pairs were saved, and the full table
was lost. Never apply a score filter at write time.

### Slots

| Slot | Type | Description |
|---|---|---|
| `$pair_scores` | data.frame | One row per tested gene pair. Columns described below. |
| `$method_params` | named list | `k` (z-score multiplier for fixed-evidence indicator, default 1.64), `weight_cap` (max pseudobulk sample count capped for weighting, default 30). |

### `$pair_scores` columns

| Column | Type | Description |
|---|---|---|
| `gene_id_A` | character | AT-ID of gene A |
| `gene_id_B` | character | AT-ID of gene B |
| `R_score` | numeric | Weighted fraction of strata with evidence. R_score = Σ(w_s · I_s) / Σ(w_s). Ranges 0–1; higher = more strata agree. |
| `z_bar` | numeric | Random-effects weighted mean Fisher z (atanh of Spearman r). |
| `tau2` | numeric | Between-stratum heterogeneity (variance component from random-effects model). |
| `pval` | numeric | Analytic weighted Poisson-binomial upper-tail p-value. Null: per-stratum null indicator probability π_s estimated from matched-permutation draws; no Monte Carlo floor. |
| `qval` | numeric | BH-FDR adjusted p-value. |
| `I_Mock` | integer (0/1) | Evidence indicator for the Mock stratum: 1 if z_s ≥ k·SE_s, else 0. One column per stratum level; column name is "I_" + stratum name. |
| `I_DC3000` | integer (0/1) | Evidence indicator for DC3000. (example; actual columns depend on stratum levels in config.) |
| `I_AvrRpt2` | integer (0/1) | Evidence indicator for AvrRpt2. |
| `I_AvrRpm1` | integer (0/1) | Evidence indicator for AvrRpm1. |
| `star` | logical | TRUE if this pair also has R_score ≥ threshold in a second dataset (cross-dataset replication annotation; added by a separate cross-dataset function). NA if cross-dataset data not available. |

**How to read R_score:** A pair with R_score = 0.75 was robustly co-expressed in
75% of the weighted stratum evidence. This is not the same as a correlation
value — it is a weighted vote count.

**How to read I_s columns:** Each I column is 1 (evidence in that condition) or 0
(no evidence). A pair that is robustly co-expressed in Mock and DC3000 but not in
the ETI conditions would have I_Mock=1, I_DC3000=1, I_AvrRpt2=0, I_AvrRpm1=0,
R_score=0.5 (assuming equal weights).

---

## ModuleInput

**One-line description:** The assembled input to the interpretation layer:
module assignments, per-module summaries, hub genes, TF members, and eigengenes.

**Produced by:** `R/interpret.R` (assembles from NetworkResult or RobustnessResult)  
**Consumed by:** `R/goi_lookup.R`, `R/interpret.R` (internally for enrichment steps)  
**Saved to disk:**

- `results/{dataset_id}/modules/module_input.rds` — full R object (RDS)
- Individual tables also saved as CSV; see below.

### Slots

| Slot | Type | Description |
|---|---|---|
| `$gene_module` | data.frame | One row per gene (all expressed genes, e.g. ~24,670). Saved as `gene_module.csv`. |
| `$module_meta` | data.frame | One row per module. Saved as `module_meta.csv`. |
| `$module_hier` | data.frame | Sub-module to top-module nesting map. Saved as `module_hier.csv`. Empty data.frame for Louvain sets (no hierarchical structure). |
| `$hub_genes` | data.frame | Top hub genes per module. Saved as `hub_genes.csv`. |
| `$module_tfs` | data.frame | TF members per module. Saved as `module_tfs.csv`. |
| `$eigengenes` | matrix (samples × modules) | Module eigengenes. Saved as `eigengenes.csv`. |
| `$method` | character | Module construction method: `"wgcna_p1"` \| `"louvain"` \| user-defined. Set by the calling script, not by the core `build_*` functions. |
| `$graph` | character | Graph label: `"large"` \| `"small"` \| user-defined. Corresponds to the R_score threshold used to filter the edge set. |
| `$r_score_threshold` | numeric | The R_score minimum applied to `rob$pair_scores` before module construction (e.g. `0.5` for the large graph, `0.6` for the small graph on 4-condition pathogen data). |

### `$gene_module` columns

| Column | Type | Description |
|---|---|---|
| `gene_id` | character | AT-ID (Araport11) |
| `top_module` | integer | Coarse-level module assignment. Module 0 = unassigned (WGCNA "grey" module). |
| `sub_module` | integer | Fine-level sub-module assignment within the top module. |
| `kME` | numeric | Module membership: Pearson correlation between this gene's expression profile and the module eigengene. Higher = more central to the module. Range −1 to 1. |

### `$module_meta` columns

| Column | Type | Description |
|---|---|---|
| `module_id` | integer | Module identifier (matches `top_module` in gene_module). |
| `n_genes` | integer | Number of genes assigned to this module. |
| `label` | character | Plain-language biological label derived from GO terms, curated anchors, and TF content. |
| `top_organ_or_condition` | character | The organ (dev atlas) or treatment condition (pathogen) where this module has its highest eigengene activity. |
| `delta_treatment` | numeric | Difference in eigengene score between the treatment condition with highest activity and Mock. Positive = induced; negative = repressed. |
| `go_top` | character | Top GO Biological Process term by enrichment in this module. |
| `zsummary` | numeric | WGCNA modulePreservation Zsummary, or fallback mean-intramodular-|cor| z-score. Higher = better preserved across contexts. Zsummary > 10 is considered high preservation. |
| `preservation_method` | character | "wgcna" if WGCNA modulePreservation was used; "fallback_meancor" if it timed out (as occurred on the pathogen side in the CZL run on large matrices). |

### `$module_hier` columns

| Column | Type | Description |
|---|---|---|
| `sub_module` | integer | Fine-level module ID |
| `top_module` | integer | Parent coarse module ID |

### `$hub_genes` columns

| Column | Type | Description |
|---|---|---|
| `module_id` | integer | |
| `gene_id` | character | AT-ID |
| `gene_symbol` | character | Gene symbol (from annotation mapping) |
| `kME` | numeric | Module membership score |
| `hub_rank` | integer | Rank within module (1 = highest kME = strongest hub) |

### `$module_tfs` columns

Source: `Athaliana_motifs_metadata.tsv` (673 TFs; `motif_id` column = AT-ID).
Use this file, not PlantTFDB (PlantTFDB auto-download failed in the CZL run).

| Column | Type | Description |
|---|---|---|
| `module_id` | integer | |
| `gene_id` | character | AT-ID (= `motif_id` in TF metadata file) |
| `gene_symbol` | character | |
| `tf_family` | character | TF family classification from the motif metadata |

---

## GOI Lookup Table

**One-line description:** The primary collaborator-facing deliverable. For each
gene in the GOI list, it reports which module the gene is in, how central it is,
and what it is co-expressed with.

**Produced by:** `R/goi_lookup.R`  
**Consumed by:** Collaborators (wet-lab biologists, CZL, external labs)  
**Saved to disk:** `results/{dataset_id}/goi/{list_name}_lookup.csv` (CSV)

One row per GOI gene. If a GOI gene is not in the network (not expressed or
filtered), it appears with all numeric columns as NA and a note in `notes`.

### Columns

| Column | Type | Description |
|---|---|---|
| `gene_id` | character | AT-ID (Araport11) |
| `gene_symbol` | character | Gene symbol |
| `module` | integer | Top-level module assignment. 0 = unassigned. |
| `kME` | numeric | Module membership score (how central this gene is to its module). Higher = more central. |
| `hub_flag` | logical | TRUE if this gene is a hub gene for its module (i.e. appears in `$hub_genes`). |
| `zsummary` | numeric | Preservation score for this gene's module. |
| `preservation_method` | character | "wgcna" or "fallback_meancor" — how `zsummary` was computed. |
| `top_N_coexpressed_partners` | character | Semicolon-separated list of the top N co-expressed partners with weights, e.g. "AT1G23456 (0.87); AT2G34567 (0.81)". N is configurable (default 10). Partners are from the same network's edge table, ordered by weight descending. |
| `notes` | character | Any flags: "not in network", "symbol not mapped to AT-ID", etc. Empty string if no issues. |

**How a collaborator uses this table:**
A wet-lab collaborator brings a list of genes they care about (receptors,
defence regulators, etc.). They open `goi_lookup.csv` and for each gene they
can see: which co-expression module it falls in (`module`), how central it is
(`kME`, `hub_flag`), whether the module is preserved across conditions
(`zsummary`), and the top genes it is co-expressed with
(`top_N_coexpressed_partners`). This directly generates hypotheses for follow-up
experiments — e.g. "this receptor is in the same module as three known immune
regulators, co-expressed primarily under DC3000 infection, with a preserved
sub-module; good candidate for epistasis assay."

---

## File Layout Summary

```
results/
└── {dataset_id}/
    ├── network/
    │   └── {stratum_id}/
    │       ├── edge_table.csv          # NetworkResult$edge_table
    │       └── params.json             # NetworkResult$params
    ├── robustness/
    │   ├── pair_scores_full.csv        # RobustnessResult$pair_scores (ALL pairs)
    │   └── robustness_result.rds       # full RobustnessResult object
    ├── modules/
    │   ├── module_input.rds            # full ModuleInput object
    │   ├── gene_module.csv             # ModuleInput$gene_module
    │   ├── module_meta.csv             # ModuleInput$module_meta
    │   ├── module_hier.csv             # ModuleInput$module_hier
    │   ├── hub_genes.csv               # ModuleInput$hub_genes
    │   ├── module_tfs.csv              # ModuleInput$module_tfs
    │   └── eigengenes.csv              # ModuleInput$eigengenes
    └── goi/
        └── {list_name}_lookup.csv      # GOI lookup table
```

For the SingleCellGGM rerun (pathogen data), network outputs go to:

```
output_per_condition/
├── Mock/edge_table.csv
├── DC3000/edge_table.csv
├── AvrRpt2/edge_table.csv
└── AvrRpm1/edge_table.csv
```
