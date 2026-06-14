# Architecture — Extended Gene Coexpression Analysis Pipeline

## 1. Data Flow

```
config.yaml  +  plugins/
      │
      ▼
┌─────────────────────────────────────────┐
│  Input Adapter  (R/adapter_seurat.R)    │
│  Seurat object → InputBundle            │
│  • log-normalised counts matrix         │
│  • cell metadata + stratum spec         │
│  • gene meta (AT-IDs mapped)            │
│  • dataset_id                           │
└─────────────────┬───────────────────────┘
                  │  InputBundle
        ┌─────────┴──────────┐
        ▼                    ▼
┌───────────────┐  ┌──────────────────────┐
│  Pseudobulk   │  │  SingleCellGGM       │
│  estimate_    │  │  estimate_           │
│  pseudobulk.R │  │  singlecellggm.R     │
└───────┬───────┘  └──────────┬───────────┘
        │                     │
        └──────────┬──────────┘
                   │  list of NetworkResult
                   │  (one per stratum level)
                   ▼
        ┌──────────────────────┐
        │  Robustness layer    │  (optional)
        │  R/robustness.R      │
        │  R_score, z_bar, τ²  │
        │  cross-dataset star  │
        └──────────┬───────────┘
                   │  RobustnessResult → assembled as ModuleInput
                   ▼
        ┌──────────────────────┐
        │  Interpretation      │
        │  R/interpret.R       │
        │  WGCNA modules,      │
        │  preservation, hubs, │
        │  GO/TF enrichment    │
        └──────────┬───────────┘
                   │
                   ▼
        ┌──────────────────────┐
        │  GOI Lookup          │
        │  R/goi_lookup.R      │
        │  per-gene table for  │
        │  collaborators       │
        └──────────────────────┘
                   │
                   ▼
              results/
              (CSVs + RDS artefacts)
```

---

## 2. Input Adapter Contract — `R/adapter_seurat.R`

The adapter is the **only** Seurat-aware file in the package. All downstream
code receives an `InputBundle` and must never import Seurat.

### InputBundle (R list)

```r
InputBundle <- list(
  counts       = <matrix: genes × cells>,   # log-normalised; full gene universe
  cell_meta    = <data.frame>,              # see columns below
  gene_meta    = <data.frame>,              # see columns below
  stratum_spec = <named list>,              # see fields below
  dataset_id   = <character(1)>            # e.g. "pathogen_multiome"
)
```

**`$counts`** — log-normalised RNA count matrix (genes as rows, cells as
columns). Genes that pass the expression filter (default: detected in ≥ 10
cells within the requested stratum subset). Row names are AT-IDs (Araport11).
Source: `RNA` assay, `counts` slot → library-size normalise → log1p.

**`$cell_meta`** columns:

| column | type | description |
|---|---|---|
| `cell_id` | character | cell barcode (unique) |
| `stratum_var` | character | value of the stratum variable for this cell (e.g. "Mock") |
| `group_var` | character | pseudobulk grouping (e.g. cluster × sample label) |
| ... | | any additional columns from the Seurat metadata are passed through |

**`$gene_meta`** columns:

| column | type | description |
|---|---|---|
| `gene_id` | character | AT-ID (Araport11), e.g. "AT1G01010" |
| `gene_symbol` | character | gene symbol; NA if no mapping found |

**`$stratum_spec`** fields:

| field | type | description |
|---|---|---|
| `variable` | character | column name in `cell_meta` used for stratification, e.g. "condition" |
| `levels` | character vector | ordered levels to iterate over, e.g. `c("Mock","DC3000","AvrRpt2","AvrRpm1")` |

### Adapter responsibilities

1. Extract log-normalised count matrix from the correct Seurat slot (RNA / counts → normalise fresh; never from SCTransform residuals).
2. Map gene symbols to AT-IDs via the annotation object from `_config_multiome.R` or a standalone mapping table supplied via config. Genes with no AT-ID mapping are dropped with a warning that counts unmapped genes.
3. Subset cells to the requested stratum level when `stratum_level` is passed (needed for per-condition GGM); returns all cells when `stratum_level = NULL`.
4. Enforce gene expression filter: keep genes detected in ≥ `min_cells` cells (default 10) within the target cell set.
5. Return `InputBundle`. Never return a Seurat object downstream of the adapter.

---

## 3. Estimation Mode Interface

Both modes accept an `InputBundle` and return a `NetworkResult` (or a list of
`NetworkResult`, one per stratum level).

### NetworkResult (R list)

```r
NetworkResult <- list(
  edge_table  = <data.frame>,    # gene_id_A, gene_id_B, weight
  gene_ids    = <character>,     # AT-IDs of all genes in the network
  stratum_id  = <character(1)>,  # which stratum, e.g. "Mock"
  mode        = <character(1)>,  # "pseudobulk" | "singlecellggm"
  params      = <named list>,    # all parameters used (pcor_cutoff, n_iter, …)
  timestamp   = <POSIXct>
)
```

**`$edge_table`** columns: `gene_id_A` (AT-ID), `gene_id_B` (AT-ID), `weight`
(pcor for GGM; Spearman r for pseudobulk). No gene symbols in this table.

### Pseudobulk mode — `R/estimate_pseudobulk.R`

- Aggregates cells to pseudobulk by `group_var` within each stratum level.
- Computes per-stratum Spearman via rank-transform-then-Pearson.
- Returns a **list** of `NetworkResult`, one per stratum level.

### SingleCellGGM mode — `R/estimate_singlecellggm.R`

- Runs separately per stratum level (one network per condition; never pooled).
- Adapter subsets cells to the stratum level before the GGM is invoked.
- **n_iter = 100** (configurable), subsample = 2,000 genes per iteration drawn
  uniformly at random from the full gene universe.
- **Aggregation = minimum |pcor| across iterations** (the defining feature of
  the method; must not use mean, median, or last-iteration).
- pcor_cutoff = 0.02, min_cells = 10 (both configurable; Phase 0 validated 0.02).
- Gene IDs are AT-IDs throughout: the adapter maps symbols before passing the
  matrix; no symbol leakage into `edge_table`.
- Returns a **list** of `NetworkResult`, one per stratum level.

---

## 4. Robustness Layer — `R/robustness.R` (optional)

**Input:** `list` of `NetworkResult` (one per stratum level from the same dataset)  
**Output:** `RobustnessResult`

```r
RobustnessResult <- list(
  pair_scores  = <data.frame>,   # see columns below — ALL tested pairs
  method_params = <named list>   # k, weight_cap, …
)
```

**`$pair_scores`** columns:

| column | type | description |
|---|---|---|
| `gene_id_A` | character | AT-ID |
| `gene_id_B` | character | AT-ID |
| `R_score` | numeric | weighted fraction of strata with evidence: Σ w_s·I_s / Σ w_s |
| `z_bar` | numeric | random-effects weighted mean Fisher z |
| `tau2` | numeric | between-stratum heterogeneity |
| `pval` | numeric | analytic weighted Poisson-binomial upper tail |
| `qval` | numeric | BH-FDR adjusted p |
| `I_<stratum>` | integer (0/1) | per-stratum evidence indicator, one column per level |

**CRITICAL:** Save **all** tested pairs, not just filtered edges. The prior CZL
run saved only R ≥ 0.7 edges and the full pair table was lost. This must not
recur. Apply filters only at interpretation time, never at write time.

**Method parameters** (via `method_params`):

| param | default | description |
|---|---|---|
| `k` | 1.64 | z-score multiplier for fixed-evidence indicator I_s |
| `weight_cap` | 30 | max pseudobulk sample count for weight w_s = sqrt(min(n_s,30)−3) |

### Cross-dataset replication ("star" annotation)

A separate function takes two `RobustnessResult` objects (e.g. dev + pathogen)
and adds a `star` column (both R_score ≥ threshold). This is an annotation, not
a filter — do not gate the output on it.

---

## 5. Shared Output Schema — Interpretation Layer

The interpretation layer (`R/interpret.R` + `R/goi_lookup.R`) always receives a
`ModuleInput`, assembled from the output of either estimation mode.

### ModuleInput (R list)

```r
ModuleInput <- list(
  gene_module  = <data.frame>,   # gene_id, top_module, sub_module, kME
  module_meta  = <data.frame>,   # per-module summary — see columns
  module_hier  = <data.frame>,   # sub_module → top_module nesting
  hub_genes    = <data.frame>,   # top hub genes per module
  module_tfs   = <data.frame>,   # TF members per module
  eigengenes   = <matrix>        # samples × modules
)
```

**`$gene_module`** — all expressed genes (e.g. ~24,670 in the CZL run):

| column | type | description |
|---|---|---|
| `gene_id` | character | AT-ID |
| `top_module` | integer | coarse-level module ID |
| `sub_module` | integer | fine-level sub-module ID |
| `kME` | numeric | module membership (correlation of gene with module eigengene) |

**`$module_meta`** — one row per module:

| column | type | description |
|---|---|---|
| `module_id` | integer | |
| `n_genes` | integer | number of genes assigned |
| `label` | character | plain-language label derived from interpretation |
| `top_organ_or_condition` | character | highest eigengene context |
| `delta_treatment` | numeric | Δ eigengene (infected − Mock) |
| `go_top` | character | top GO BP term |
| `zsummary` | numeric | WGCNA modulePreservation Zsummary (or fallback meancor z-score) |
| `preservation_method` | character | "wgcna" \| "fallback_meancor" |

Note: `preservation_method = "fallback_meancor"` is expected when
`modulePreservation` times out on large matrices (as occurred on the pathogen
side in the CZL run).

**`$module_hier`** — sub-module to top-module nesting:

| column | type | description |
|---|---|---|
| `sub_module` | integer | fine-level module ID |
| `top_module` | integer | parent coarse module ID |

**`$hub_genes`** — top hub genes:

| column | type | description |
|---|---|---|
| `module_id` | integer | |
| `gene_id` | character | AT-ID |
| `gene_symbol` | character | |
| `kME` | numeric | |
| `hub_rank` | integer | rank within module (1 = highest kME) |

**`$module_tfs`** — TF members per module. Source: `Athaliana_motifs_metadata.tsv`
(673 TFs; use this file, NOT PlantTFDB — PlantTFDB auto-download failed in the
CZL run):

| column | type | description |
|---|---|---|
| `module_id` | integer | |
| `gene_id` | character | AT-ID (= `motif_id` column in TF file) |
| `gene_symbol` | character | |
| `tf_family` | character | TF family from motif metadata |

**`$eigengenes`** — matrix: samples (rows) × modules (columns).

### GOI Lookup — `R/goi_lookup.R`

**Input:** `ModuleInput` + a GOI list (from config `plugins:` block)  
**Output:** one-row-per-gene data.frame, the **primary collaborator-facing deliverable**.

| column | type | description |
|---|---|---|
| `gene_id` | character | AT-ID |
| `gene_symbol` | character | |
| `module` | integer | top_module assignment |
| `kME` | numeric | |
| `hub_flag` | logical | TRUE if gene is a hub gene for its module |
| `zsummary` | numeric | module preservation score |
| `preservation_method` | character | "wgcna" \| "fallback_meancor" |
| `top_N_coexpressed_partners` | character | semicolon-joined list of top N partners with weights, e.g. "AT1G23456 (0.87); AT2G34567 (0.81)" |

Design note: lab members and collaborators arrive with a list of genes they care
about; `goi_lookup.R` tells them which module each gene sits in, how central it
is, what it is co-expressed with, and in what condition context. This is the
primary deliverable. Design this function first-class, not as an afterthought.

---

## 6. Plugin Boundary

Everything dataset- or lab-specific lives in `plugins/` or is referenced via
the config YAML. Core functions must not import anything from `plugins/` directly.
Plugins are injected via config or function arguments only.

Plugin examples:

| Plugin | Description |
|---|---|
| GOI gene lists | CZL receptors/ligands, lab focus genes, collaborator lists |
| Curated anchor sets | NLR, PTI_receptor, vascular, abscission, … |
| `Athaliana_motifs_metadata.tsv` | 673 TFs; `motif_id` column = AT-ID |
| Symbol→AT-ID mapping table | used by adapter + GGM output reconciliation |

---

## 7. File Responsibilities

| File | Role |
|---|---|
| `R/adapter_seurat.R` | Seurat → InputBundle; only Seurat-aware file |
| `R/estimate_pseudobulk.R` | InputBundle → list of NetworkResult (pseudobulk) |
| `R/estimate_singlecellggm.R` | InputBundle → list of NetworkResult (GGM) |
| `R/robustness.R` | list of NetworkResult → RobustnessResult |
| `R/interpret.R` | RobustnessResult / NetworkResult → ModuleInput |
| `R/goi_lookup.R` | ModuleInput + GOI list → GOI table (collaborator deliverable) |
| `config/example_config.yaml` | Reference config with all settable parameters |
| `plugins/` | Dataset- and lab-specific assets; not imported by core |
| `inst/scripts/run_pipeline.R` | Top-level driver; reads config; calls adapter → estimation → robustness → interpret → goi_lookup |

---

## Module construction methods (under evaluation)

For GGM mode, the interpretation layer currently supports multiple module
construction paths:

- **WGCNA** at explicit soft powers (1, 4, 6, 8) — `build_wgcna_modules()`
  with `soft_power` override. Auto power selection on GGM output converges
  to power=1 with poor scale-free fit (FLAG-11).
- **Louvain** — `igraph::cluster_louvain()` on the R_score-filtered weighted
  graph (weight = abs(tanh(z_bar))).
- **Leiden** — `igraph::cluster_leiden()` with objective_function='modularity'
  on the same graph.

A comprehensive benchmarking layer (5 thresholds × 6 methods = 30 cells)
evaluates structure-only metrics (modularity, grey rate, module sizes, Gini,
cross-method ARI/NMI) without biological pre-judgment.
See `inst/scripts/benchmark_modules_pathogen.R` and
`results/pathogen_multiome/method_benchmark/`.

This benchmarking step may become a permanent, general-purpose
`benchmark_module_methods()` pipeline component — see BENCHMARK_REPORT.md
design note section for the proposed API.

---

## Official module sets — pathogen multiome

The benchmark resolved to a **dual-method strategy** on **two graphs**, producing
four official module sets as standard output:

| Set | Graph | Method | R_score threshold |
|---|---|---|---|
| `large_wgcna` | large | WGCNA power=1 | ≥ 0.5 |
| `large_louvain` | large | Louvain | ≥ 0.5 |
| `small_wgcna` | small | WGCNA power=1 | ≥ 0.6 |
| `small_louvain` | small | Louvain | ≥ 0.6 |

### Why two methods

- **WGCNA power=1** — conservative, hierarchical (top + sub module structure),
  high inter-module separation, suitable for identifying compact co-regulatory
  cores. `sub_module` is meaningful only for WGCNA sets.
- **Louvain** — comprehensive, maximises modularity across the full graph,
  tends to produce more and larger modules with higher assigned-gene coverage.
  `sub_module = NA` throughout (no hierarchical structure); `module_hier` is an
  empty data.frame for Louvain sets.

Both methods run on both graphs as standard output. Which set to use for a given
downstream analysis is left to the user/analyst.

### Why two graphs (R_score thresholds)

On 4-condition pathogen data, R_score is **discrete** with only 5 possible
values: 0, 0.25, 0.5, 0.75, 1.0 (corresponding to 0/1/2/3/4 conditions meeting
the evidence criterion). Thresholds 0.3 and 0.4 therefore resolve to the same
graph as 0.5; thresholds 0.6 and 0.7 resolve to the same graph as 0.6.
Consequently there are exactly **two distinct graphs**:

- **Large graph** (R_score ≥ 0.5): ~62,863 edges, ~10,358 genes — pairs
  robust in at least 2/4 conditions.
- **Small graph** (R_score ≥ 0.6): ~15,384 edges, ~3,441 genes — pairs
  robust in at least 3/4 conditions.

On future datasets with more strata (e.g. 10 organs), R_score is continuous and
threshold selection becomes a continuous dial rather than this discrete step
function. The threshold logic in `build_wgcna_modules()` and the Louvain helper
accepts any numeric value; only the effective graph count changes.

### ModuleInput method identifier fields

`ModuleInput` now carries three additional metadata fields (set at construction
time by the official-module script; not produced by the core `build_*` functions
themselves):

| Field | Type | Values |
|---|---|---|
| `$method` | character | `"wgcna_p1"` \| `"louvain"` \| user-defined |
| `$graph` | character | `"large"` \| `"small"` \| user-defined label |
| `$r_score_threshold` | numeric | e.g. `0.5`, `0.6` |

See also `docs/OUTPUT_SCHEMA.md` for the updated ModuleInput slot table.

---

## Condition-pattern characterization

`characterize_condition_pattern()` in `R/robustness.R` extends the robustness
layer with a flat, comprehensive description of each pair's cross-condition
activity. It operates on a `RobustnessResult` plus the per-condition
`NetworkResult` list and returns one row per pair with two output tracks:

### Discrete pattern (threshold-based)

The four `I_<condition>` indicator columns already produced by
`compute_robustness()` are concatenated into a 4-character bit string in
`condition_order` (default: Mock / DC3000 / AvrRpt2 / AvrRpm1). All 16 possible
patterns are assigned a `pattern_label`:

| Pattern | Label |
|---|---|
| `1111` | `constitutive_all` |
| `0111` | `pan_pathogen` |
| `0011` | `ETI_shared` |
| `0001` | `single_AvrRpm1` |
| `0010` | `single_AvrRpt2` |
| `0100` | `single_DC3000` |
| `1000` | `single_Mock` |
| `0000` | `none` |
| all others | `mixed_<pattern>` |

**These labels are mechanical names for bit patterns.** They carry no
claim about regulatory mechanism, epistasis, or causal structure. "ETI_shared"
means bits 3 and 4 are both set; it does not claim the genes are involved in
ETI. Use the labels as grouping handles, not biological assertions.

### Continuous profile (threshold-free)

For each condition, the raw edge weight (pcor) is looked up from the
per-condition `edge_table`. Pairs absent from a condition's network receive
weight 0. Derived columns: `w_max`, `w_min`, `w_range`, `w_mean`, and a
`specificity_index`:

```
specificity_index = (w_max - w_mean_of_others) / (w_max + ε)
  where w_mean_of_others = (Σ w - w_max) / (S - 1),  ε = 1e-6
```

Ranges from 0 (uniform across all conditions) toward 1 (signal concentrated in
a single condition).

### Usage

The full pair table (`pair_condition_patterns.csv`, all 1.4M pairs, FLAG-03
compliant) supports both flat filtering by `pattern_label` and continuous
clustering on the `w_<condition>` matrix. The companion
`module_condition_patterns.csv` (one row per module per set) summarises the
distribution of patterns within each module (`dominant_pattern`,
`frac_<label>`, mean `w_<condition>`, `module_specificity_index`).

This is the recommended way to find condition-specific modules: filter or sort
`all_modules_condition_patterns.csv` by `dominant_pattern` or
`frac_single_<condition>` rather than by R_score, which structurally
de-prioritises condition-specific pairs.

### Generalisation

The function is not specific to the pathogen 4-condition design. Any
`RobustnessResult` with multiple strata will produce a meaningful pattern
characterisation. For S strata, the pattern string has S bits; named shortcuts
(`pan_pathogen`, `ETI_shared`) apply only to the default 4-condition order.

---

## Estimation mode and rare cell populations

A structural comparison of the two estimation modes on the same gene can guide
mode selection for rare-population genes.

### Empirical finding (AT1G08860 / BON3, pathogen multiome)

- **Cell-level SingleCellGGM (pooled):** AT1G08860 was placed in the grey
  (unassigned) module with near-zero kME and approximately 11 weak,
  condition-private partners (each appearing in a single condition only).
- **Subcluster-level pseudobulk:** AT1G08860 was placed in an assigned module
  with kME ≈ 0.60 and approximately 1,032 partners.

### Structural explanation

GGM conditions on ~18k genes using all cells simultaneously. For a gene whose
expression is elevated only in a rare cell population, the signal is diluted by
the large number of cells where the gene is near zero. The partial-correlation
estimate reflects the marginal cell-level structure, which is dominated by
abundant cell types.

Pseudobulk at the subcluster level uses subclusters as observations. A
subcluster in which the gene is highly expressed becomes one data point that
drives the marginal Spearman correlation with co-expressed partners. Rare
populations therefore contribute in proportion to the number of subclusters they
span, not in proportion to their cell count — recovering co-expression that
cell-level GGM cannot detect.

### Guidance for mode selection

- **Pseudobulk (subcluster grouping)** is the recommended first mode for
  recovering co-expression of rare-cell-population genes. The number of
  subclusters observed determines statistical power, not the raw cell count
  for the rare type.
- **Cell-level SingleCellGGM** is better suited for dense, broadly-expressed
  network structure across the full gene universe and for removing confounding
  by indirect regulation (partial vs marginal correlation). Its effective
  resolution degrades for genes expressed in a small fraction of cells.

The two modes are complementary, not competing. For a comprehensive analysis
of a dataset with both broad and rare populations, running both and comparing
the partner lists and module assignments for genes of interest is informative.
