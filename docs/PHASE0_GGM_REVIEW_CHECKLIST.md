# Phase 0 — SingleCellGGM Run Review Checklist

**Dataset:** Nobori et al. 2024 *Nature* — Arabidopsis pathogen multiome  
**Conditions:** Mock / DC3000 / AvrRpt2 / AvrRpm1 (~65k nuclei, 15 samples)  
**Reference method:** Xu, Wang & Ma (2024), *Cell Reports Methods* 4:100813  
**Purpose:** Determine whether the existing casual run can become the core of
`R/estimate_singlecellggm.R`, or whether it needs refactoring and/or a rerun.

Each item is answered: **Pass / Fail / NA** + one-line evidence.

---

## Axis 1 — Parameter Validity

### pcor cutoff

- [ ] **1.1** The pcor threshold used is ≥ 0.03 (paper default), OR a different value is
  explicitly stated with justification.  
  _Evidence:_ `threshold = ___`

- [ ] **1.2** The cutoff is applied to the final conservative |pcor| (minimum across
  iterations), not to per-iteration pcor values.  
  _Evidence:_

### min-cells threshold

- [ ] **1.3** The minimum number of cells in which an edge must appear is ≥ 10 (paper
  default), OR a different value is stated with justification.  
  _Evidence:_ `min_cells = ___`

- [ ] **1.4** "Appearing in a cell" is defined correctly: the edge (non-zero partial
  correlation between gene pair) is counted per cell, not per subsample or per iteration.  
  _Evidence:_

### Subsampling iterations and aggregation rule

- [ ] **1.5** The number of subsampling iterations is recorded (paper uses enough for
  convergence; typically ≥ 100 is expected).  
  _Evidence:_ `n_iterations = ___`

- [ ] **1.6** The aggregation rule is **minimum |pcor| across iterations** (conservative;
  the defining feature of SingleCellGGM). A different aggregation (mean, median, max,
  last-iteration) would mean the run does not implement the method as published.  
  _Evidence:_ aggregation used = `___`

- [ ] **1.7** The 2,000-gene subsample per iteration is drawn **uniformly at random from
  the full gene universe** (not stratified, not pre-filtered to a gene set of interest).  
  _Evidence:_

### Gene universe

- [ ] **1.8** The number of genes in the input matrix is recorded.  
  _Evidence:_ `n_genes = ___`

- [ ] **1.9** Genes were not pre-filtered to a biologically specific subset (e.g. LR genes,
  receptor genes) before running — such filtering would bias the GEPs and invalidate
  generalization.  
  _Evidence:_

- [ ] **1.10** Low-expressed genes were filtered by an explicit, recorded criterion
  (e.g. expressed in ≥ N cells, mean count ≥ X). The criterion is recoverable.  
  _Evidence:_ filter used = `___`

### Input matrix normalization

- [ ] **1.11** The normalization/transform applied to the cell-level matrix is recorded
  (log-normalized? SCTransform? raw counts? other?).  
  _Evidence:_ transform = `___`

- [ ] **1.12** The choice is consistent with SingleCellGGM's assumptions. The method was
  benchmarked on log-normalized data; SCTransform residuals or raw counts are deviations
  that need justification.  
  _Evidence:_

### Cell selection and filtering

- [ ] **1.13** The number of cells actually fed to the model is recorded (all ~65k, or a
  subset?).  
  _Evidence:_ `n_cells_used = ___`

- [ ] **1.14** Doublet removal status is recorded (were called doublets excluded?).  
  _Evidence:_

- [ ] **1.15** Any additional cell quality filters (min UMI, min genes, % mito) are recorded
  and their thresholds are the same as those used for the published analysis, or deviations
  are noted.  
  _Evidence:_

### Condition handling

- [ ] **1.16** It is explicitly documented whether conditions (Mock/DC3000/AvrRpt2/AvrRpm1)
  were **pooled into one run** or **run separately per condition**.  
  _Evidence:_ mode = `pooled | per-condition | ___`

- [ ] **1.17** The implication is recorded: a pooled run produces a single network that
  conflates condition effects (treatment response and co-expression are entangled); a
  per-condition run produces four networks amenable to cross-condition comparison, which is
  consistent with the pipeline's per-stratum architecture.  
  _Evidence:_

- [ ] **1.18** If pooled: it is acknowledged that condition identity is a major source of
  variation and that the pooled GEPs may reflect condition clusters rather than co-regulatory
  programs.  
  _Evidence:_

---

## Axis 2 — Output Structure

### Artifacts present

- [ ] **2.1** An edge list (gene1, gene2, pcor) or equivalent sparse pcor matrix exists
  on disk and its path is recorded.  
  _Evidence:_ path = `___`

- [ ] **2.2** GEP / module assignments (gene → module ID) exist, OR a clustering step was
  not run and only the raw network is present (note which).  
  _Evidence:_

- [ ] **2.3** Gene-level metadata (kME equivalent, module membership scores, hub flags)
  exist, or it is noted they are absent and must be computed downstream.  
  _Evidence:_

- [ ] **2.4** Per-cell scores (cell loadings on each GEP, if computed) exist, or it is
  noted they are absent.  
  _Evidence:_

### Schema compatibility with pipeline output contract

- [ ] **2.5** The edge list / pcor matrix can be coerced to the pipeline's shared network
  format (genes × genes sparse matrix, or long-format gene1/gene2/weight data.frame)
  without information loss.  
  _Evidence:_

- [ ] **2.6** Module assignments can be coerced to the pipeline's module contract
  (named integer vector: gene_id → module_id) as consumed by `interpret.R`.  
  _Evidence:_

- [ ] **2.7** Items explicitly missing from the current output that the pipeline requires
  are listed (e.g. kME, preservation stats, GO enrichment — these are expected to be absent
  at this stage and computed later).  
  _Evidence:_ missing = `___`

### Gene ID format

- [ ] **2.8** Gene identifiers in the output are **AT-IDs (Araport11 / TAIR10 AGI format,
  e.g. AT1G01010)**. Mixed formats (symbol + AT-ID, or only symbol) require reconciliation
  before the GOI lookup resource can be built.  
  _Evidence:_ ID format = `___`

- [ ] **2.9** Gene IDs are compatible with the dev atlas gene universe (GSE226097) for
  future cross-dataset work. If the pathogen run used a different annotation build or
  filtered to a different gene universe, the overlap must be checked.  
  _Evidence:_

### GOI lookup recoverability

- [ ] **2.10** Per-gene module membership is recoverable from the existing output, so that
  `build_goi_table()` can be populated for any arbitrary GOI list without re-running the
  GGM.  
  _Evidence:_

---

## Axis 3 — Reproducibility

### Randomness and determinism

- [ ] **3.1** A random seed is set and recorded before the subsampling loop. Without a
  seed, the GEPs are not exactly reproducible.  
  _Evidence:_ `set.seed(___)`

- [ ] **3.2** The GEPs are stable enough across re-runs to serve as a reliable reference
  (if seed is absent or convergence was not checked, this is a known gap to flag).  
  _Evidence:_

### Software environment

- [ ] **3.3** The SingleCellGGM package version (or commit hash if installed from GitHub)
  is recorded.  
  _Evidence:_ `SingleCellGGM v___`

- [ ] **3.4** R version and key dependency versions (Seurat, Matrix, etc.) are recorded
  (sessionInfo() or renv.lock or equivalent).  
  _Evidence:_

### Input recoverability

- [ ] **3.5** The exact input Seurat object (or count matrix) used is identified by an
  unambiguous path or checksum — not "the object I had loaded in my session at the time."  
  _Evidence:_ path/checksum = `___`

- [ ] **3.6** The cell selection (which of the ~65k nuclei were used) is reproducible from
  the script alone, without relying on session state or intermediate objects that may no
  longer exist.  
  _Evidence:_

### Runtime and rerunability

- [ ] **3.7** Approximate wall-clock runtime for the run is noted (this is the heavy mode).  
  _Evidence:_ runtime ≈ `___`

- [ ] **3.8** The run is rerunnable from the script as-is on the lab machine without
  manual intervention (correct paths, no missing intermediate files, no session-state
  dependencies).  
  _Evidence:_

---

## Axis 4 — Generalizability (can this become the pipeline core?)

### Hardcoding

- [ ] **4.1** The script contains no hardcoded paths that are specific to the pathogen
  dataset and would silently produce wrong results on a different input (e.g. hardcoded
  gene lists, cell barcodes, sample names, condition labels).  
  _Evidence:_

- [ ] **4.2** The script contains no hardcoded gene counts or cell counts used as
  parameters (e.g. subsetting to exactly 65k cells).  
  _Evidence:_

- [ ] **4.3** The script does not assume a specific number of conditions or specific
  condition labels (Mock/DC3000/etc.) — these must be driven by the config's stratum spec.  
  _Evidence:_

### Input adapter compatibility

- [ ] **4.4** The script accepts a normalized count matrix + cell metadata as inputs,
  OR it can be refactored to accept them without restructuring the core algorithm.
  If it reads directly from a Seurat object in a way that is not isolated, note the
  refactor needed.  
  _Evidence:_

- [ ] **4.5** The script would work on the dev atlas (GSE226097) or a future dataset
  without modification beyond supplying a different input matrix and stratum spec.  
  _Evidence:_

### Refactor scope

- [ ] **4.6** The refactoring required to wrap this run behind `R/estimate_singlecellggm.R`
  is estimated and categorized:
  - [ ] **Minimal** — parameterize paths + seed + thresholds; no algorithmic changes
  - [ ] **Moderate** — condition-handling logic needs restructuring; input adapter needed
  - [ ] **Substantial** — aggregation rule is wrong or normalization is incompatible;
    rerun required before wrapping  
  _Evidence:_

### Known method limitation (record explicitly — not a failure criterion)

- [ ] **4.7** It is explicitly recorded in the review notes that SingleCellGGM — like all
  co-expression methods — **cannot recover paracrine ligand–receptor pairs**. Single-cell
  resolution makes paracrine pairs *harder* to detect, not easier (ligand and receptor in
  different cells → pcor ≈ 0 or negative at cell level). This run is **not** evaluated on
  LR-pair recovery. The pipeline's value is co-regulatory module discovery and context
  interpretation.  
  _Evidence:_ acknowledged = `yes / no`

---

## VERDICT

Fill in after completing all items above.

| Category | Count |
|---|---|
| Pass | |
| Fail | |
| NA | |

### Decision rubric

**Use as-is (rare):** All Axis 1 parameters match the paper, aggregation rule is correct,
conditions handled per-stratum, seed set, env captured, input recoverable, no hardcoding.
Wrap directly into `R/estimate_singlecellggm.R` with minimal parameterization.

**Refactor, no rerun needed:** Parameters are valid and aggregation rule is correct, but
the script is hardcoded or session-dependent. Refactor for generalizability and input-adapter
compatibility; existing output artifacts are trustworthy.

**Rerun required:** Aggregation rule deviates from minimum-across-iterations (Fail 1.6),
OR normalization is incompatible (Fail 1.12), OR conditions were pooled when per-condition
networks are needed for the pipeline's per-stratum architecture (Fail 1.16/1.18) and pooling
was not an intentional design choice. Rerun with corrected parameters before wrapping.

**Verdict (fill in):**

> ___ [ use as-is / refactor, no rerun / rerun required ]
>
> Key blockers (if any): ___
>
> Recommended next action: ___
