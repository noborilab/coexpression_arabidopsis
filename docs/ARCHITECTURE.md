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
