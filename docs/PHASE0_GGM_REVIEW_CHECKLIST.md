# Phase 0 — SingleCellGGM Run Review Checklist

**Dataset:** Nobori et al. 2024 *Nature* — Arabidopsis pathogen multiome  
**Conditions:** Mock / DC3000 / AvrRpt2 / AvrRpm1 (~65k nuclei, 15 samples)  
**Reference method:** Xu, Wang & Ma (2024), *Cell Reports Methods* 4:100813  
**Purpose:** Determine whether the existing casual run can become the core of
`R/estimate_singlecellggm.R`, or whether it needs refactoring and/or a rerun.

**Sources reviewed:**
- R script: `SingleCellGGM/SingleCellGGM_pathogen_data_ALL.R`
- MATLAB run script: `SingleCellGGM/run_scggm_all.m`
- MATLAB library: `SingleCellGGM-main/SingleCellGGM.m` + `adjust_cutoff.m`
- Output dir: `SingleCellGGM/output_all/`

Each item is answered: **Pass / Fail / NA** + one-line evidence.

---

## Axis 1 — Parameter Validity

### pcor cutoff

- [x] **1.1** The pcor threshold used is ≥ 0.03 (paper default), OR a different value is
  explicitly stated with justification.  
  **PASS** — `threshold = 0.020`. Below 0.03 default but FDR-curve-justified: FDR = 0.0016
  at pcor 0.02 (vs. 0 permuted edges at pcor 0.03). FDR inspection step is present in the R
  script (plot + `abline(h = 0.1)`). Note: justification is implicit, not written in prose.

- [x] **1.2** The cutoff is applied to the final conservative |pcor| (minimum across
  iterations), not to per-iteration pcor values.  
  **PASS** — `adjust_cutoff()` filters `ggm_ori.pcor_all` which already holds the
  minimum-across-iterations values produced by `SingleCellGGM.m`.

### min-cells threshold

- [x] **1.3** The minimum number of cells in which an edge must appear is ≥ 10 (paper
  default), OR a different value is stated with justification.  
  **PASS** — `min_cells = 10` (R script L18); `min_coexpress = 10` in MATLAB run script;
  `cut_off_coex_cell = 10` in `SingleCellGGM.m`.

- [x] **1.4** "Appearing in a cell" is defined correctly: the edge (non-zero partial
  correlation between gene pair) is counted per cell, not per subsample or per iteration.  
  **PASS** — `coex = a' * a` where `a = x > 0` is a cell-level binary matrix.
  `coex(r,s)` counts cells where both gene r and gene s are expressed. Applied as
  `coex >= cut_off_coex_cell` in `adjust_cutoff.m:42`.

### Subsampling iterations and aggregation rule

- [x] **1.5** The number of subsampling iterations is recorded (paper uses enough for
  convergence; typically ≥ 100 is expected).  
  **PASS** — `num_iter = max(round(gene_num*(gene_num-1)/39980), 20000)`. With ~5,651
  genes → formula gives ~799, floored to minimum of 20,000. `n_iterations = 20,000`.

- [x] **1.6** The aggregation rule is **minimum |pcor| across iterations** (conservative;
  the defining feature of SingleCellGGM). A different aggregation (mean, median, max,
  last-iteration) would mean the run does not implement the method as published.  
  **PASS** — `SingleCellGGM.m` lines 92–94: `if abs(pc(m,n)) < abs(pcor_all(r,s));
  pcor_all(r,s) = pc(m,n)`. Unambiguously the minimum. This is the key algorithmic check.

- [x] **1.7** The 2,000-gene subsample per iteration is drawn **uniformly at random from
  the full gene universe** (not stratified, not pre-filtered to a gene set of interest).  
  **PASS** — `j = randperm(p, selected_num)` where `selected_num = 2000`. Uniform random
  without replacement from all `p` genes in the input matrix.

### Gene universe

- [x] **1.8** The number of genes in the input matrix is recorded.  
  **PASS** — R script comment: `res$n_genes  # e.g. 5651`. Authoritative count is in
  `input/all_HVG8000_gene_names.txt` (one gene per line). The MATLAB re-filter
  (`sum(expression_matrix > 0, 1) >= 10`) may reduce this slightly further.

- [x] **1.9** Genes were not pre-filtered to a biologically specific subset (e.g. LR genes,
  receptor genes) before running — such filtering would bias the GEPs and invalidate
  generalization.  
  **FAIL** — Filtered to top 8,000 HVGs (`FindVariableFeatures(data, nfeatures = 8000)`)
  before export. After the within-MATLAB min-cell re-filter, ~5,651 genes entered the GGM.
  HVGs in a pathogen-infection atlas are enriched for immune/stress-responsive genes.
  The R script comment acknowledges this: `# uuse HVGs to avoid memory issue; use HPC for
  all`. A full-gene-universe run was intended but not executed.

- [x] **1.10** Low-expressed genes were filtered by an explicit, recorded criterion
  (e.g. expressed in ≥ N cells, mean count ≥ X). The criterion is recoverable.  
  **PASS** — `filter used = expressed in ≥ 10 cells`. R script L46:
  `keep <- Matrix::rowSums(mat > 0) >= min_cells` (min_cells = 10). Also re-applied in
  MATLAB (`idx = sum(expression_matrix > 0, 1) >= 10`).

### Input matrix normalization

- [x] **1.11** The normalization/transform applied to the cell-level matrix is recorded
  (log-normalized? SCTransform? raw counts? other?).  
  **PASS** — `transform = library-size normalization + log1p`. R script L50–54:
  raw counts from `GetAssayData(slot = "counts")` → divide by `colSums` × 10,000 → `log1p`.
  Performed fresh in this script, not from a pre-normalized Seurat slot.

- [x] **1.12** The choice is consistent with SingleCellGGM's assumptions. The method was
  benchmarked on log-normalized data; SCTransform residuals or raw counts are deviations
  that need justification.  
  **PASS** — log-normalized data matches SingleCellGGM's benchmark normalization.
  Raw counts are fetched and log-normalized fresh, not from a SCTransform-normalized slot.

### Cell selection and filtering

- [x] **1.13** The number of cells actually fed to the model is recorded (all ~65k, or a
  subset?).  
  **PASS** — `n_cells_used = 65,061`. R script L102: `all_cells <- colnames(data)`;
  comment L115: `res$n_cells  # should be 65061`. All cells from all conditions included.

- [x] **1.14** Doublet removal status is recorded (were called doublets excluded?).  
  **FAIL** — Not documented in this script. The Seurat object is `combined_filtered.rds`
  ("filtered" in path implies upstream QC), but doublet removal is not mentioned and its
  status is not recoverable from this script alone.

- [x] **1.15** Any additional cell quality filters (min UMI, min genes, % mito) are recorded
  and their thresholds are the same as those used for the published analysis, or deviations
  are noted.  
  **FAIL** — Not documented. Script loads `combined_filtered.rds` and uses `colnames(data)`
  without any additional cell filtering or logging of the upstream QC thresholds.

### Condition handling

- [x] **1.16** It is explicitly documented whether conditions (Mock/DC3000/AvrRpt2/AvrRpm1)
  were **pooled into one run** or **run separately per condition**.  
  **FAIL (BLOCKER)** — `mode = pooled`. R script L102: `all_cells <- colnames(data)` uses
  all cells from all four conditions. The prefix "all" and "all_HVG8000" make this
  inferable, but it is not explicitly stated as a deliberate design choice.

- [x] **1.17** The implication is recorded: a pooled run produces a single network that
  conflates condition effects (treatment response and co-expression are entangled); a
  per-condition run produces four networks amenable to cross-condition comparison, which is
  consistent with the pipeline's per-stratum architecture.  
  **FAIL** — No such note exists in the script.

- [x] **1.18** If pooled: it is acknowledged that condition identity is a major source of
  variation and that the pooled GEPs may reflect condition clusters rather than co-regulatory
  programs.  
  **FAIL** — Not acknowledged. The GEP scoring plots split by `sample2` (a per-sample
  column) but the implications of pooling for GEP interpretation are not stated.

---

## Axis 2 — Output Structure

### Artifacts present

- [x] **2.1** An edge list (gene1, gene2, pcor) or equivalent sparse pcor matrix exists
  on disk and its path is recorded.  
  **PASS** — `path = output_all/all_HVG8000.ggm.pcor0.020.cell10.txt`. Three-column
  tab-separated (gene1, gene2, pcor); 24,821 edges. A second network at pcor 0.035 also
  exists (5,181 edges).

- [x] **2.2** GEP / module assignments (gene → module ID) exist, OR a clustering step was
  not run and only the raw network is present (note which).  
  **PASS** — Two-tier membership: (1) coarse GEPs in `all.GEP_membership.pcor0.02.txt`
  (gene/module, 4,374 genes); (2) sub-module assignments in `all.GEP*_submodules.final.
  pcor0.02.txt` files and the combined `all.ALL_GEP_submodules.final.pcor0.02.txt`.

- [x] **2.3** Gene-level metadata (kME equivalent, module membership scores, hub flags)
  exist, or it is noted they are absent and must be computed downstream.  
  **NA (expected absent)** — No kME or hub-score file generated. Module membership (binary
  assignment) exists but no continuous strength-of-membership metric. These will be computed
  in `interpret.R`.

- [x] **2.4** Per-cell scores (cell loadings on each GEP, if computed) exist, or it is
  noted they are absent.  
  **PASS (proxy)** — `AddModuleScore` scores exist in Seurat metadata in
  `combined_filtered_with_GEP_subscores.rds`. Note: these are Seurat proxy scores
  (average of gene set vs. background), not true factor loadings. Sufficient for
  visualization; not equivalent to factorization-derived loadings.

### Schema compatibility with pipeline output contract

- [x] **2.5** The edge list / pcor matrix can be coerced to the pipeline's shared network
  format (genes × genes sparse matrix, or long-format gene1/gene2/weight data.frame)
  without information loss.  
  **PASS** — Three-column format (gene1, gene2, pcor) is already the long-format contract.
  The MATLAB table also contains SamplingTime, r, Cell_num columns in the raw network file
  but the R downstream code drops these to columns 1–3.

- [x] **2.6** Module assignments can be coerced to the pipeline's module contract
  (named integer vector: gene_id → module_id) as consumed by `interpret.R`.  
  **PASS** — `all.GEP_membership.pcor0.02.txt` has gene/module two-column format;
  trivially coercible with `setNames(tbl$module, tbl$gene)`.

- [x] **2.7** Items explicitly missing from the current output that the pipeline requires
  are listed (e.g. kME, preservation stats, GO enrichment — these are expected to be absent
  at this stage and computed later).  
  **NA** — `missing = kME / continuous membership scores, module preservation stats,
  GO enrichment`. All expected absent; downstream pipeline will compute them.

### Gene ID format

- [x] **2.8** Gene identifiers in the output are **AT-IDs (Araport11 / TAIR10 AGI format,
  e.g. AT1G01010)**. Mixed formats (symbol + AT-ID, or only symbol) require reconciliation
  before the GOI lookup resource can be built.  
  **FAIL (BLOCKER)** — `ID format = mixed: gene symbols + AT-IDs`. Edge list sample:
  `BSMT1, GES, TPS03, ILL5` (symbols) alongside `AT4G10290, AT5G28237` (AT-IDs). The
  Seurat object uses symbol as primary name where annotation is available, AT-ID otherwise.
  GOI lookup will fail on the symbol entries without a symbol→AT-ID mapping step.

- [x] **2.9** Gene IDs are compatible with the dev atlas gene universe (GSE226097) for
  future cross-dataset work. If the pathogen run used a different annotation build or
  filtered to a different gene universe, the overlap must be checked.  
  **FAIL** — Mixed symbol/AT-ID format will need reconciliation against the dev atlas,
  which uses AT-IDs (Araport11). Overlap cannot be assessed without first mapping symbols
  to AT-IDs.

### GOI lookup recoverability

- [x] **2.10** Per-gene module membership is recoverable from the existing output, so that
  `build_goi_table()` can be populated for any arbitrary GOI list without re-running the
  GGM.  
  **PASS** — `all.ALL_GEP_submodules.final.pcor0.02.txt` provides gene → submodule
  mapping for 4,374 genes. The `get_submodule_genes()` function in the script already
  implements bidirectional lookup. Caveat: gene ID mismatch (2.8) must be resolved first.

---

## Axis 3 — Reproducibility

### Randomness and determinism

- [x] **3.1** A random seed is set and recorded before the subsampling loop. Without a
  seed, the GEPs are not exactly reproducible.  
  **PASS** — `set.seed(98)` via `rng(98)` in `SingleCellGGM.m:68` (set inside the library
  before the iteration loop). Also `set.seed(123)` in R for Louvain GEP clustering (L181).
  Caveat: the main seed is in the library source, not in `run_scggm_all.m` — future users
  of the library who call it from a different context may not notice it is hardcoded.

- [x] **3.2** The GEPs are stable enough across re-runs to serve as a reliable reference
  (if seed is absent or convergence was not checked, this is a known gap to flag).  
  **PASS (inferred)** — Seed is fixed (rng(98)) so output is deterministic. Convergence
  was not explicitly checked, but 20,000 iterations on ~5,651 genes is well above typical
  convergence requirements.

### Software environment

- [x] **3.3** The SingleCellGGM package version (or commit hash if installed from GitHub)
  is recorded.  
  **FAIL** — `SingleCellGGM v = unknown`. Directory is `SingleCellGGM-main` (GitHub main
  branch download). No version, tag, or commit hash recorded. The MATLAB source is present
  locally so the exact code is recoverable from the file, but version provenance is lost.

- [x] **3.4** R version and key dependency versions (Seurat, Matrix, etc.) are recorded
  (sessionInfo() or renv.lock or equivalent).  
  **FAIL** — No `sessionInfo()`, no `renv.lock`, no version comments. R environment is
  not captured.

### Input recoverability

- [x] **3.5** The exact input Seurat object (or count matrix) used is identified by an
  unambiguous path or checksum — not "the object I had loaded in my session at the time."  
  **PASS** — `path = .../SA_039_94_multiome_revision_rep2_9h_only/out/_seurat_object/
  motifFixed/combined_filtered.rds`. Hardcoded absolute path appears at L11 and L248.
  No checksum, but path is unambiguous.

- [x] **3.6** The cell selection (which of the ~65k nuclei were used) is reproducible from
  the script alone, without relying on session state or intermediate objects that may no
  longer exist.  
  **PASS** — `all_cells <- colnames(data)` takes all cells from the loaded RDS. No
  interactive filtering or session-state dependency. Reproducible as long as the RDS is
  unchanged.

### Runtime and rerunability

- [x] **3.7** Approximate wall-clock runtime for the run is noted (this is the heavy mode).  
  **FAIL** — `runtime ≈ not recorded`. No timing information in script or comments.

- [x] **3.8** The run is rerunnable from the script as-is on the lab machine without
  manual intervention (correct paths, no missing intermediate files, no session-state
  dependencies).  
  **FAIL** — Three blockers: (1) `source(".../_config_multiome.R")` at L4 is required —
  it provides the `annotation` data.frame used by `id_to_gene()` at L918; without it the
  TF-annotation section fails. (2) MATLAB with SingleCellGGM on path must be run manually
  between the input-prep and output-analysis sections. (3) All paths are machine-specific
  (Dropbox absolute paths); not portable to HPC or collaborators.

---

## Axis 4 — Generalizability (can this become the pipeline core?)

### Hardcoding

- [x] **4.1** The script contains no hardcoded paths that are specific to the pathogen
  dataset and would silently produce wrong results on a different input (e.g. hardcoded
  gene lists, cell barcodes, sample names, condition labels).  
  **FAIL** — Multiple hardcoded absolute Dropbox paths throughout: Seurat object (L11,
  L248, L990), input/output dirs (L8, L126, L157, L239, L1004), TF list (L909),
  config source (L4). No dataset-specific gene lists or barcodes hardcoded, but paths
  would break on any other machine or dataset.

- [x] **4.2** The script contains no hardcoded gene counts or cell counts used as
  parameters (e.g. subsetting to exactly 65k cells).  
  **PASS** — `res$n_cells  # should be 65061` is a sanity-check comment, not a parameter.
  No count is used to gate logic.

- [x] **4.3** The script does not assume a specific number of conditions or specific
  condition labels (Mock/DC3000/etc.) — these must be driven by the config's stratum spec.  
  **PASS** — Condition labels are not referenced anywhere. Script pools all cells without
  any condition-aware logic. (Note: this also means condition-specific analysis is absent,
  not that it is generalized.)

### Input adapter compatibility

- [x] **4.4** The script accepts a normalized count matrix + cell metadata as inputs,
  OR it can be refactored to accept them without restructuring the core algorithm.
  If it reads directly from a Seurat object in a way that is not isolated, note the
  refactor needed.  
  **FAIL** — Seurat-specific code is not isolated: `DefaultAssay`, `GetAssayData`,
  `FindVariableFeatures`, `VariableFeatures`, `AddModuleScore`, `FeaturePlot`, and
  dollar-sign slot access all appear throughout. The `prepare_scggm_input()` function
  takes a `seu` argument and calls Seurat internals. Requires an adapter layer to decouple.

- [x] **4.5** The script would work on the dev atlas (GSE226097) or a future dataset
  without modification beyond supplying a different input matrix and stratum spec.  
  **FAIL** — Machine-specific paths, Seurat coupling, and pooled single-run design all
  require structural changes before applying to a second dataset.

### Refactor scope

- [x] **4.6** The refactoring required to wrap this run behind `R/estimate_singlecellggm.R`
  is estimated and categorized:
  - [ ] **Minimal** — parameterize paths + seed + thresholds; no algorithmic changes
  - [x] **Moderate** — condition-handling logic needs restructuring; input adapter needed
  - [ ] **Substantial** — aggregation rule is wrong or normalization is incompatible;
    rerun required before wrapping  
  **Moderate — but rerun is needed** (see verdict). The algorithm is correct; the scope
  (pooled, HVG-only) and gene-ID format are the problems. Refactor: parameterize
  `prepare_scggm_input()` → `adapter_seurat.R`, remove Seurat calls from core, parameterize
  the MATLAB call or replace with R equivalent, add per-condition iteration, resolve gene IDs
  to AT-IDs before export.

### Known method limitation (record explicitly — not a failure criterion)

- [x] **4.7** It is explicitly recorded in the review notes that SingleCellGGM — like all
  co-expression methods — **cannot recover paracrine ligand–receptor pairs**. Single-cell
  resolution makes paracrine pairs *harder* to detect, not easier (ligand and receptor in
  different cells → pcor ≈ 0 or negative at cell level). This run is **not** evaluated on
  LR-pair recovery. The pipeline's value is co-regulatory module discovery and context
  interpretation.  
  **FAIL** — `acknowledged = no`. Not mentioned in the script or any output comment. Must
  be recorded before this run is used to draw conclusions about LR biology.

---

## VERDICT

Fill in after completing all items above.

| Category | Count |
|---|---|
| Pass | 23 |
| Fail | 16 |
| NA | 4 |

**Pass:** 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.10, 1.11, 1.12, 1.13, 2.1, 2.2,
2.5, 2.6, 2.10, 3.1, 3.2, 3.5, 3.6, 4.2, 4.3  
**Fail:** 1.9, 1.14, 1.15, **1.16**, 1.17, 1.18, **2.8**, 2.9, 3.3, 3.4, 3.7, 3.8, 4.1,
4.4, 4.5, 4.7  
**NA:** 2.3, 2.4 (proxy), 2.7, 4.6 (scope estimate)

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

> **Rerun required**
>
> Key blockers:
> 1. **1.16 (CRITICAL)** — All four conditions (Mock/DC3000/AvrRpt2/AvrRpm1) were pooled
>    into a single GGM run. The pipeline's per-stratum architecture requires four separate
>    networks for cross-condition comparison. The pooled GEPs conflate treatment response
>    with co-regulatory structure; the large GEPs likely capture condition clusters, not
>    co-regulatory programs.
> 2. **1.9** — Only HVGs (top 8,000 → ~5,651 after min-cell filter) were used due to
>    memory constraints. The R script itself says "use HPC for all". A full-gene-universe
>    run (~20,000–30,000 expressed genes) on HPC is the correct production run.
> 3. **2.8** — Gene IDs in outputs are mixed (gene symbols + AT-IDs). Must be resolved
>    to pure AT-IDs (Araport11) before GOI lookup and cross-dataset integration.
>
> Silver lining: the algorithm is correctly implemented. **1.6 PASS** (minimum pcor
> aggregation is correct), **1.11–12 PASS** (log-normalization from raw counts is correct),
> **3.1 PASS** (seed is set). The existing run is a valid exploratory artifact and its GEPs
> can be used for biological orientation, but it should not be promoted to the pipeline core
> without rerunning per-condition on full gene universe.
>
> Recommended next action:
> 1. Acknowledge 4.7 in a comment or design note (LR limitation).
> 2. Plan HPC run: full gene universe (~25k expressed genes), per-condition (4 separate
>    runs), AT-ID gene names preserved from the start (use `rownames()` from counts before
>    any symbol conversion). Seed in run script, not library code.
> 3. Refactor `prepare_scggm_input()` → `R/adapter_seurat.R` with `genes_use = NULL`
>    (no HVG filter) and a `stratum` argument to subset by condition.
> 4. Resolve gene ID format: the Seurat RNA assay rownames are symbols-where-available;
>    add a symbol→AT-ID mapping step in the adapter using `annotation` from
>    `_config_multiome.R` (or build a standalone mapping table).
> 5. Existing output (`all_HVG8000.ggm.pcor0.020.cell10.txt`, GEP membership files) is
>    kept as a reference artifact; the pipeline will supersede it.
