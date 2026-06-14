# Pipeline Design Flags & Decisions

Running log of issues found during Phase 0–1 that affect implementation or
interpretation. Update each flag's Status line as work progresses.

---

## FLAG-01: n_iter default
**Phase**: 1 → corrected Phase 2c
**Issue**: The casual pathogen GGM run used ~20,000 iterations (formula-driven).
Pipeline spec originally set n_iter = 100, misreading the paper: 100 is the
average number of times each GENE PAIR is sampled, not the total iteration count.
With n_iter = 100 on a real gene universe (p ~ 16k), most pairs are never sampled.
**Decision**: CORRECTED. n_iter = NULL (auto) computes round(p*(p-1)/39980) per
stratum so each pair is sampled ~100x on average. For p = 16k this gives ~6,400
iterations. The formula is implemented in estimate_singlecellggm.R.
**Status**: CORRECTED in Phase 2c. See estimate_singlecellggm.R: resolved_n_iter
block in the per-stratum loop.

---

## FLAG-02: preservation_method = "fallback_meancor"
**Phase**: 1
**Issue**: WGCNA modulePreservation timed out on the pathogen side in the CZL run
(~1,500 samples x 11k genes). The lightweight fallback (mean intramodular |cor|
z-score) was used instead. Zsummary values in module_preservation.csv may reflect
the fallback, not the full WGCNA Zsummary.
**Decision**: Document via preservation_method column ("wgcna" | "fallback_meancor").
interpret.R propagates this field; no special-casing needed.
**Status**: Schema captures it. Implement gracefully in Phase 2e (interpret.R).

---

## FLAG-03: Full R_score pair table not saved in CZL run
**Phase**: 0
**Issue**: CZL run saved only R >= 0.7 filtered edges (14,497 rows). The raw
per-pair R_score table for all ~25k tested pairs was lost — exists only inside
run-checkpoint RDS files.
**Decision**: robustness.R MUST save pair_scores_full.csv (all tested pairs) before
any filtering. Enforced in RobustnessResult schema.
**Status**: Schema enforces it. Verify in Phase 2d (robustness.R).

---

## FLAG-04: Mixed gene IDs in casual GGM output
**Phase**: 0
**Issue**: Edge list from the casual run contains gene symbols (BSMT1, ILL5, TPS03)
mixed with AT-IDs (AT4G10290). Breaks GOI lookup against Araport11 AT-IDs.
**Decision**: adapter_seurat.R maps symbols -> AT-IDs at input time using symbol_map
argument. No symbols enter the edge_table. Unmapped genes logged with warning count.
**Status**: Implemented in Phase 2a (adapter_seurat.R).

---

## FLAG-05: Conditions pooled in casual GGM run — BLOCKER
**Phase**: 0
**Issue**: All 65,061 nuclei from Mock + DC3000 + AvrRpt2 + AvrRpm1 were fed into
a single GGM run. GEPs likely encode condition identity, not co-regulatory programs.
**Decision**: Rerun per condition (4 separate NetworkResults). Plan in docs/RERUN_PLAN.md.
**Status**: Casual run artifacts NOT usable as pipeline core. Rerun required
(parallel to Phase 2 implementation).

---

## FLAG-06: HVG filter in casual GGM run — BLOCKER
**Phase**: 0
**Issue**: FindVariableFeatures(nfeatures = 8000) reduced input to ~5,651 genes.
Script comment "# use HPC for all" confirms full-universe run was never done.
**Decision**: Rerun with all expressed genes (~25k), no HVG filter.
Expression filter: detected in >= 10 cells per condition subset.
**Status**: Rerun required. Parameters in docs/RERUN_PLAN.md.

---

## FLAG-07: GOI lookup partners — long format companion needed
**Phase**: 1
**Issue**: top_N_coexpressed_partners is a semicolon-joined string in the lookup CSV.
Easy to read in Excel but hard to use programmatically.
**Decision**: Add goi_partners_long.csv as a companion output from goi_lookup.R.
Schema: gene_id, partner_id, partner_symbol, weight, rank
**Status**: Open. Implement in Phase 2e (goi_lookup.R).

---

## FLAG-08: Dev atlas pseudobulk expression matrix not found on disk
**Phase**: 1
**Issue**: The 428 sub-clusters x 24,670 genes expression matrix that was the WGCNA
input is not in any checked directory. chk_v2_shared.rds (85 MB) likely contains
eigengenes but not necessarily the full matrix. Cell-level Seurat objects for dev
atlas are archived off-machine.
**Decision**: Not blocking Phase 2a-2c. Revisit when pseudobulk mode is wired to
the dev atlas.
**Status**: Open.

---

## FLAG-09: min_cells / coex_cutoff semantics
**Phase**: 2b → corrected Phase 2c
**Issue**: The Phase 2b spec and the original implementation conflated two distinct
quantities under the name `min_cells`:
- (a) number of *iterations* in which a gene pair was sampled (sampling count)
- (b) number of *cells* in which both genes are co-detected (count > 0 in same cell)
The original code filtered on (a), using the iteration count as a proxy for (b).
These are unrelated: a pair can be sampled in many iterations but co-detected in
very few cells (or vice versa).
**Decision**: CORRECTED. The parameter is renamed `coex_cutoff` and now correctly
filters on the co-detection cell count, computed from a p×p co-detection matrix
`coex = tcrossprod(counts > 0)` built once per stratum before the iteration loop.
The sampling count (samp > 0) is kept as a separate filter to exclude never-sampled
pairs. See estimate_singlecellggm.R: coex matrix construction and the keep_edge
filter in the post-loop section.
**Status**: CORRECTED in Phase 2c.

---

## FLAG-10: Sign handling — positive-only vs |pcor| filter
**Phase**: 2c
**Issue**: The original implementation retained edges where |pcor| >= pcor_cutoff,
keeping both positive and negative partial correlations. The SingleCellGGM paper
(Xu et al. 2024) retains only positive partial correlations (pcor >= cutoff), as
co-expression networks represent co-activation rather than mutual inhibition.
**Decision**: Pipeline default `keep_negative = FALSE` matches the paper (positive
pcor only). A toggle `keep_negative = TRUE` is exposed for users who want signed
networks or wish to study negative partial correlations. See estimate_singlecellggm.R:
the keep_edge filter in the post-loop section.
**Status**: Implemented in Phase 2c.

---

## FLAG-11: WGCNA soft-thresholding on GGM partial-correlation networks
**Phase**: benchmark (2026-06)
**Issue**: WGCNA soft-thresholding was designed for dense correlation matrices.
On GGM output (already sparse; partial correlations), auto soft-power selection
chose power=1 and scale-free R² fit was poor (max R²≈0.67). The benchmark
(inst/scripts/benchmark_modules_pathogen.R) compared WGCNA at powers 1/4/6/8
and graph-clustering methods (Louvain, Leiden) across five R_score thresholds
(0.3–0.7) on the pathogen multiome GGM robustness results. Evaluation was
structure-only: modularity, grey rate, module size, cross-method ARI.
**Decision**: RESOLVED — pipeline adopts dual-method strategy: WGCNA power=1
(conservative, high inter-module separation, hierarchical sub-module structure)
+ Louvain (comprehensive, high modularity, no hierarchy). Both methods run on
both graphs (large: R_score ≥ 0.5; small: R_score ≥ 0.6) as standard output.
The choice of which result set to use for a given downstream analysis is left to
the user/analyst. See inst/scripts/run_official_modules_pathogen.R for the
implementation and docs/ARCHITECTURE.md for the full rationale.
**Status**: RESOLVED (2026-06). Four official module sets produced.

---

## FLAG-12: R_score consistency vs condition-specificity
**Phase**: post-benchmark (2026-06)
**Issue**: The R_score weighted-consistency statistic structurally favours
all-condition pairs (I_s = 1 in all strata → R_score = 1.0) and buries
condition-specific pairs (e.g. a pair present only in AvrRpm1 gets R_score ≈
0.25). But condition-specific modules are often of primary biological interest.
Filtering on R_score therefore suppresses exactly the subset users most want to
examine when studying condition-specific responses.
**Decision**: Added `characterize_condition_pattern()` in `R/robustness.R`,
producing:
- **Discrete 4-bit pattern** per pair (16 possible patterns; bit order = condition_order):
  named labels for the 8 most interpretable patterns ("constitutive_all",
  "pan_pathogen", "ETI_shared", "single_<condition>", "none") and
  "mixed_<pattern>" for the remaining 8. Labels are mechanical — they name bit
  patterns, not biological mechanisms.
- **Continuous per-condition weights** (w_Mock, w_DC3000, w_AvrRpt2, w_AvrRpm1)
  drawn directly from the per-condition edge tables, plus w_max, w_min, w_range,
  w_mean, and a specificity_index = (w_max − w_mean_of_others) / (w_max + ε).
All 1,413,505 pairs are characterised and saved (FLAG-03 compliant).
R_score is retained as one descriptor among several, not the primary ranking.
`inst/scripts/compute_condition_patterns_pathogen.R` generates the pair table;
`inst/scripts/module_condition_patterns.R` profiles each official module set.
**Status**: IMPLEMENTED (2026-06-14).

---

## FLAG-13: GGM misses inducible/cell-specific transcription factors; pseudobulk captures them
**Phase**: post-benchmark (2026-06)
**Issue**: Cell-level SingleCellGGM misses inducible or cell-type-restricted transcription
factors, as shown at family scale on the WRKY TF family (70 genes): GGM placed only 14/70
in a confident module (4 more at near-zero/negative kME); subcluster pseudobulk captured
47 additional WRKYs that GGM missed, many at high kME (e.g. WRKY8 kME 0.92 / 2858 partners;
WRKY75 kME 0.89; WRKY35 kME 0.94). This reproduces the single-gene BON3 finding (FLAG/architecture
note on rare cell populations) at family scale. Structural cause: GGM conditions over all cells
and all genes, diluting pulse-induced or cell-type-restricted signal; pseudobulk uses subclusters
as observations.
**Decision**: GUIDANCE: for inducible or cell-type-restricted regulators, use pseudobulk
(subcluster grouping). The two modes are complementary. Both modes are standard pipeline output.
**Status**: documented; both modes are standard pipeline output.

---

## FLAG-14: Observation-point design as a first-class variable
**Phase**: Stage 0-2 (2026-06-14)
**Issue**: Standard clustering optimises cell classification; co-expression needs
observation points that spread along covariation axes and denoise without
destroying within-group signal. Treating the observation-point design as a
tunable variable — rather than a fixed implementation detail — reveals a
fundamental degree of freedom in pseudobulk co-expression analysis. The
design choice (granularity, aggregation method, normalization) has a large
effect on which genes are "visible" and how stable the resulting network is.
**Decision**: Implement multiple in-house observation-point generators
(`obs_cluster`, `obs_subcluster`, `obs_metacell_knn`, `obs_stratified`,
`obs_axis_bin`) with a standardised `ObsPointSet` interface, and a prior-free
evaluation harness (`coexpr_eval.R`) that scores designs without reference to
GO terms or known gene sets. The evaluation Pareto front (stability vs richness)
guides design selection empirically. Comparison methods (hdWGCNA, CS-CORE,
SuperCell, SEACells) are deferred; the interface is method-agnostic so they can
be slotted in later. TF annotation remains post-hoc interpretation only — never
edge estimation (see GRN boundary note in docs/ARCHITECTURE.md).
**Status**: Stage 0-2 implemented (2026-06-14). Runner script:
`inst/scripts/obs_design_sweep_pathogen.R`.
