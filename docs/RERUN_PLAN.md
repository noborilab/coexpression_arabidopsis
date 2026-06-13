# SingleCellGGM Rerun Plan — Pathogen Multiome (per-condition, full gene universe)

**Status:** Data task. Can run in parallel with Phase 2 implementation.  
**This is NOT a coding task.** The algorithm is already implemented in MATLAB
(`SingleCellGGM-main/`). The task is to run it four times with corrected scope.

**Motivation (from Phase 0 review):**
- Blocker 1.16: the existing run pooled all four conditions into one GGM.
  The pipeline requires per-condition networks for cross-condition comparison.
- Blocker 1.9: only HVGs (~5,651 genes) were used. The full expressed gene
  universe (~25k genes) is required.
- Blocker 2.8: gene IDs in outputs are mixed (symbols + AT-IDs). Outputs must
  use AT-IDs throughout.

---

## 1. Input

**Seurat object path:**
```
/Users/jep23kod/Nobori Lab (TSL) Dropbox/Tatsuya NOBORI/SALK_clowd/Projects/
SA_PTI_ETI_single_cell/SA_039_94_multiome_revision_rep2_9h_only/out/
_seurat_object/motifFixed/combined_filtered.rds
```

**Which slot to use for counts:**  
RNA assay, `counts` slot (raw integer counts). Do not use SCTransform or a
pre-normalised slot. Normalise fresh: divide by library size × 10,000, then
log1p. This matches the paper benchmark normalization and the Phase 0 run.

**Which cells:**  
All 65,061 nuclei present after existing upstream QC (doublets already removed
in the `combined_filtered.rds` object; no additional cell filtering needed).
The QC thresholds used to produce this object are documented in the published
analysis; do not re-apply them here.

---

## 2. Condition Split

Run the GGM **four times**, once per condition:

| Run | Condition | Metadata column | Value to subset |
|---|---|---|---|
| 1 | Mock | (confirm column name from Seurat metadata) | "Mock" |
| 2 | DC3000 | same | "DC3000" |
| 3 | AvrRpt2 | same | "AvrRpt2" |
| 4 | AvrRpm1 | same | "AvrRpm1" |

Confirm the exact column name and level spellings by inspecting the Seurat
object metadata (`colnames(seu@meta.data)`, `table(seu$condition)`) before
running.

**Do NOT pool conditions.** The pooled run in `output_all/` is kept as a
reference artifact but will not be used by the pipeline.

---

## 3. Gene Universe

Use **all expressed genes** — no HVG filter.

Expression filter per condition subset:
- Keep genes detected (raw count > 0) in ≥ **10 cells** within that condition's
  cell subset.
- Apply this filter independently per condition (a gene absent in Mock but
  present in DC3000 will pass the DC3000 filter only).
- Expected: ~20,000–25,000 genes per condition subset after filtering.

The R script comment from the Phase 0 run already flagged this: `# use HVGs to
avoid memory issue; use HPC for all`. This is that run.

---

## 4. Gene ID Handling

**Use AT-IDs from the start.**

1. After loading the Seurat object and subsetting to one condition, check
   `rownames(GetAssayData(seu_sub, assay="RNA", slot="counts"))`. These are
   mixed: some are gene symbols, some are AT-IDs (as found in the Phase 0 review).
2. Load the annotation object from `_config_multiome.R` (or export it as a
   standalone TSV: two columns, `gene_symbol` and `gene_id` / AT-ID). This is
   the same annotation used in the original analysis.
3. Map all rownames: symbol → AT-ID. For genes that already have AT-ID format
   (AT[1-5MC]G[0-9]{5}), keep as-is. For symbols, look up in the annotation
   table.
4. Log the count of unmapped genes (symbols with no AT-ID match) before
   filtering. Drop unmapped genes and add the count to the run log.
5. Rename matrix rows to AT-IDs before any further processing.
6. The MATLAB input matrix row names must be AT-IDs. The output edge list column
   names (`gene1`, `gene2`) will inherit these AT-IDs. Do not re-map after the
   GGM — the IDs must be AT-IDs from entry to exit.

---

## 5. Parameters

Carry forward from Phase 0 validation. These are confirmed correct.

| Parameter | Value | Notes |
|---|---|---|
| `n_iter` | 100 | Default for R-based re-implementation. The Phase 0 MATLAB run used a formula yielding ~20,000 iterations, which was overspecified for convergence; 100 is the paper's intent and sufficient. |
| `subsample` | 2,000 | Genes sampled per iteration (paper default). |
| `aggregation` | minimum \|pcor\| | Conservative; the defining feature. Do not use mean/median/last. |
| `pcor_cutoff` | 0.02 | Phase 0 validated: FDR = 0.0016 at this threshold (vs. 0 permuted edges at 0.03). |
| `min_cells` | 10 | Minimum cells where both genes are co-expressed (paper default). |
| `seed` | 98 | Must be set in the run script (not buried in the library source as in the Phase 0 run). |

---

## 6. Output

Save results matching the `NetworkResult` schema (see `docs/OUTPUT_SCHEMA.md`).

**Output directory:** Use a NEW directory, distinct from `output_all/` (the Phase 0 pooled run). Do not overwrite existing output.

**Suggested path pattern:**
```
output_per_condition/
├── Mock/
│   ├── edge_table.csv     # columns: gene_id_A, gene_id_B, weight
│   └── params.json        # all run parameters as JSON
├── DC3000/
│   ├── edge_table.csv
│   └── params.json
├── AvrRpt2/
│   ├── edge_table.csv
│   └── params.json
└── AvrRpm1/
    ├── edge_table.csv
    └── params.json
```

**`edge_table.csv` columns:**

| Column | Content |
|---|---|
| `gene_id_A` | AT-ID of gene A |
| `gene_id_B` | AT-ID of gene B |
| `weight` | minimum \|pcor\| across iterations (the GGM output) |

**`params.json` content (example):**
```json
{
  "stratum_id": "Mock",
  "mode": "singlecellggm",
  "n_cells": 16241,
  "n_genes": 23814,
  "n_unmapped_genes_dropped": 12,
  "n_iter": 100,
  "subsample": 2000,
  "aggregation": "min_abs_pcor",
  "pcor_cutoff": 0.02,
  "min_cells": 10,
  "seed": 98,
  "timestamp": "2026-06-XX ..."
}
```

---

## 7. Runtime

**Unknown until first condition completes.**

Background:
- Phase 0 pooled run: ~5,651 genes, 65,061 cells, 20,000 iterations. Runtime not
  recorded (Fail 3.7 in Phase 0 review).
- This rerun: ~20,000–25,000 genes per subset, ~15,000–18,000 cells per condition
  subset, 100 iterations. The gene count is 4–5× higher; iteration count is ~200×
  lower; cell count per run is ~3–4× lower.
- Net effect on runtime is uncertain without profiling.

**Recommendation: run Mock first as a timing probe.**  
Record wall-clock time for Mock before queuing the other three conditions. If
runtime is acceptable on the MacBook Pro, run the remaining three sequentially.
If it exceeds ~4 hours per condition, move to HPC (the R script comment from
Phase 0 already anticipated this need: `# use HPC for all`).

Flag in the run log: "Runtime for this run: X minutes on [machine]."

---

## Checklist Before Running

- [ ] Confirm condition column name and level spellings in Seurat metadata
- [ ] Load and inspect `_config_multiome.R` to locate the `annotation` object
- [ ] Export or verify symbol→AT-ID mapping table (two columns: `gene_symbol`, `gene_id`)
- [ ] Verify `output_per_condition/` directory does not exist or is empty
- [ ] Confirm `SingleCellGGM-main/` MATLAB library is on path
- [ ] Set seed = 98 in the run script (not in the MATLAB library)
- [ ] Log unmapped gene count before each run
- [ ] Time the Mock run before queueing remaining conditions
