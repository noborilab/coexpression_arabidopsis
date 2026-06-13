# Architecture

## Pipeline overview

```
per-context network estimation → cross-context preservation/robustness → modules + interpretation
```

> **TODO (Phase 1):** Expand this into a full architecture document once the scaffold is reviewed
> and Phase 0 (review of existing SingleCellGGM run) is complete.

---

## Input adapter boundary

The core pipeline operates on a single abstraction:

> **(normalized counts matrix) + (cell metadata data.frame) + (stratum spec)**

`Seurat` is today's input format, handled exclusively by `R/adapter_seurat.R`.
No other file in the package imports or calls Seurat.

Future adapters (AnnData via `reticulate`, raw-count re-normalization from GEO,
direct count matrices) implement the same four-field contract and are drop-in
replacements without touching any core logic:

```
adapter_*(path, ...) → list(
  counts      = <genes × cells matrix>,
  meta        = <cell metadata data.frame>,
  stratum_spec = <character vector of context column names>,
  dataset_id  = <short string>
)
```

Long-term goal: integrate ALL published single-cell datasets, which will require
re-normalizing from raw reads. The adapter boundary ensures this never forces
rewrites of core estimation or interpretation code.

---

## Estimation modes (both first-class)

| Mode | File | When to use |
|---|---|---|
| Pseudobulk | `R/estimate_pseudobulk.R` | Default; light compute; requires pseudobulk structure in dataset |
| SingleCellGGM | `R/estimate_singlecellggm.R` | Cell-level; heavier compute; cell-level objects required |

Both modes produce the same downstream interface (gene × gene edge table or
matrix) consumed by the robustness and interpretation layers.

**Phase 0 gate:** the existing casual SingleCellGGM run on the pathogen
multiome data (Nobori 2024) must be reviewed before `estimate_network_singlecellggm()`
is implemented. Parameters to audit: pcor cutoff, min-cells threshold,
subsampling iterations, reproducibility, output format.

---

## Robustness layer (optional)

`R/robustness.R` — computes cross-stratum R_score and cross-dataset replication.

- Enabled/disabled via `config$robustness$enabled`.
- Skipped entirely when the dataset has too few independent strata to be meaningful.
- Does not gate downstream modules; robustness scores are annotations, not filters.

---

## Output / interpretation layer (shared)

`R/interpret.R` and `R/goi_lookup.R` — shared by both estimation modes.

- Module construction (WGCNA signed network, soft-power, merge threshold)
- Module hierarchy: coarse + nested sub-modules (granularity sweep)
- Cross-context preservation (WGCNA `modulePreservation` Zsummary, or
  lightweight intramodular-|cor| z-score fallback for large matrices)
- Hub genes (kME), TF intersection, GO BP enrichment, curated-set fold-enrichment
- GOI lookup table: per-gene module, kME, hub flag, preservation

---

## Plugins (swappable domain resources)

`plugins/` — domain-specific gene lists and curated sets referenced by path
in the config. **Never hardcoded** in package source.

Examples:
- `plugins/Athaliana_motifs_metadata.tsv` — TF list (673 TFs, column `motif_id`)
- `plugins/curated_anchors.rds` — curated gene-set anchors for fold-enrichment

---

## Data flow (planned)

```
config.yaml
    │
    ▼
[Input Adapter]          adapter_seurat.R  (or future adapter_*)
    │  counts + meta + stratum_spec
    ▼
[Estimation]             estimate_pseudobulk.R  OR  estimate_singlecellggm.R
    │  per-stratum networks
    ▼
[Robustness] (optional)  robustness.R
    │  edge-level R_score + replication annotations
    ▼
[Interpretation]         interpret.R  +  goi_lookup.R
    │  modules, preservation, hubs, GO/TF/anchor enrichment, GOI table
    ▼
results/
```
