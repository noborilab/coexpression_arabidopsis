# Session Handoff — 2026-06-14/15

Covers work done in the session ending ~02:00 on 2026-06-15.
Picks up from FLAG-12 (condition-pattern layer) and documents through FLAG-14 (observation-point design sweep).

---

## 1. Decisions made this session

### Normalization: zscore_gene + Spearman (FLAG-14)

Empirically settled on the pathogen multiome data (11010 genes, 298 subclusters).
Full table from `results/pathogen_multiome/obs_design/normalization_decision.csv`:

| norm_method   | cor_type | depth_leak | splithalf_cor | splithalf_jaccard | eff_rank | frac_visible |
|---|---|---|---|---|---|---|
| none          | spearman | 0.806      | 0.895         | 0.606             | 36.2     | 0.996        |
| none          | pearson  | 0.806      | 0.831         | 0.527             | 36.2     | 0.996        |
| cp10k_log     | spearman | 0.652      | 0.807         | 0.591             | 133.6    | 1.000        |
| cp10k_log     | pearson  | 0.652      | 0.733         | 0.582             | 133.6    | 1.000        |
| log_only      | spearman | 0.811      | 0.895         | 0.605             | 63.8     | 0.996        |
| log_only      | pearson  | 0.811      | 0.848         | 0.579             | 63.8     | 0.996        |
| **zscore_gene** | **spearman** | **0.107** | **0.895** | 0.605         | **159.9** | **1.000**  |
| zscore_gene   | pearson  | 0.107      | 0.832         | 0.530             | 159.9    | 1.000        |

**Selection rule**: lowest depth_leakage_rho among options whose splithalf_mat_cor is ≥ 90% of the best.
**Winner: zscore_gene + Spearman**.

Why this matters:
- **depth_leakage 0.107 vs 0.806–0.811** for raw/log-only: z-scoring per gene removes almost all depth/abundance confound, so hubs are not just highly expressed genes.
- **eff_rank 159.9 vs 36.2** for raw: z-scoring resolves 4.4× more independent covariation axes because it eliminates the single dominant "all genes up/down together" axis that depth drives, exposing the fine-grained cell-state variation underneath.
- **splithalf_cor 0.895** matches the best of all other methods: high reproducibility is preserved despite more axes being resolved.
- Spearman > Pearson for reproducibility across all normalization methods (~0.06 gap consistently).

**What this means for obs-point design**: pseudobulk co-expression is now treated as a first-class tunable variable (FLAG-14). The choice of how to aggregate cells and how to normalize the aggregates determines what structure is recoverable. Normalization is an empirically settled question for this dataset; it should be re-run on any new dataset.

---

## 2. Stage 2 granularity sweep results

**Settings**: zscore_gene + Spearman, run via `inst/scripts/obs_design_sweep_pathogen.R`.
**Note**: Metacell designs (t200/t100/t50/t25) were terminated due to runtime (PCA on 65k cells per split-half rep → ~6 hours remaining for 4 designs). Results below cover completed designs only.

### Full table (completed designs)

| design           | n_pts | eff_rank | pred_r2  | splithalf | note |
|---|---|---|---|---|---|
| **subcluster**   | **298** | **159.9** | **0.249** | **0.895** | Pareto-dominant |
| cluster_res0.10  | 34    | 23.5     | −0.155   | 0.537     | fallback bug (see below) |
| cluster_res0.25  | 34    | 23.5     | −0.155   | 0.537     | fallback bug |
| cluster_res0.50  | 34    | 23.5     | −0.155   | 0.537     | fallback bug |
| cluster_res1.00  | 30    | 18.2     | +0.256   | 0.754     | fallback bug (different column hit) |
| cluster_res2.00  | 34    | 23.5     | −0.155   | 0.537     | fallback bug |
| cluster_res4.00  | 34    | 23.5     | −0.155   | 0.537     | fallback bug |
| metacell_t200    | —     | —        | —        | —         | terminated (runtime) |
| metacell_t100    | —     | —        | —        | —         | terminated (runtime) |
| metacell_t50     | —     | —        | —        | —         | terminated (runtime) |
| metacell_t25     | —     | —        | —        | —         | terminated (runtime) |

### Pareto front (stability vs richness)

**Among completed designs, `subcluster` (298 pts) Pareto-dominates everything else simultaneously on both axes**:
- Stability (splithalf): 0.895 vs 0.537–0.754 for all cluster designs
- Richness (eff_rank): 159.9 vs 18.2–23.5 for all cluster designs; pred_r2=0.249 vs ≤ 0.256

The cluster designs all collapse to the same 34 pre-computed Seurat clusters regardless of the resolution parameter (see Bug #1). The one anomaly is cluster_res1.00 which hit a different metadata column (giving 30 pts and slightly different metrics) but still dramatically worse than subcluster.

**Conclusion (provisional)**: The 298-subcluster design is the recommended observation-point layout for this dataset. Metacell designs may be competitive but cannot be confirmed without completing the sweep. Given the subcluster is based on biologically meaningful pre-computed sub-clustering and already gives eff_rank=159.9, the prior expectation is that metacell designs would be comparable rather than strictly superior — but this needs empirical confirmation.

### obs_cluster resolution fallback: critical finding

All 6 `obs_cluster` resolution sweeps returned the same 34 points because `obs_cluster` finds any column whose name matches `"cluster"` in the Seurat metadata and reuses it, ignoring the `resolution` parameter. This made the cluster granularity sweep uninformative. **Fix required before re-running** (see Bug #1 below).

---

## 3. Next steps

### Immediate (next session)

1. **Re-run pseudobulk co-expression with zscore_gene + Spearman on the 298-subcluster design.**
   This is the Pareto-optimal design from the completed sweep. Use `obs_subcluster(bundle, group_col="sub_clst_rna_20260610")` followed by `normalize_obs(obs, "zscore_gene")` and `coexpr_from_obs(obs, cor_type="spearman")`. Wire into the existing robustness layer and module construction pipeline.

2. **Fix Bug #1 (obs_cluster resolution fallback)** before re-running the cluster granularity sweep. The fix is narrow: see Known bugs below.

3. **Complete metacell sweep** after Bug #1 fix and runtime optimization. The metacell designs need a faster PCA path (e.g. irlba) or should reduce split-half reps for the sweep (3 reps instead of 5 for the first exploration). The key scientific question is whether manifold-tiling metacells add axes beyond what subclusters already resolve.

### Medium term

4. **Feature plots**: run `inst/scripts/featureplot_modules.R` once the new pseudobulk network is available from the zscore_gene + subcluster run.

5. **BON3 / WRKY post-hoc sanity**: run the post-hoc section of `obs_design_sweep_pathogen.R` on the subcluster design once the network is available. Confirm BON3 visibility and WRKY partner recovery (this was blocked in the sweep by the metacell termination; the subcluster post-hoc section was never reached in Stage 2).

6. **Dev atlas**: when cell-level Seurat objects are restored from archive (FLAG-08), run the obs-design sweep on the dev atlas data. The same normalization question should be re-settled empirically there; zscore_gene may or may not be optimal across tissue types.

7. **Comparison methods** (hdWGCNA, CS-CORE, SuperCell, SEACells): the `ObsPointSet` interface is method-agnostic. These can be slotted in as additional `design_fn` arguments to `evaluate_obs_design` without touching the evaluation harness. Add in a later session.

---

## 4. Known bugs to fix

### Bug 1 — obs_cluster ignores resolution parameter (OPEN)

**File**: `R/observation_points.R`, function `obs_cluster`  
**Symptom**: All 6 resolution values (0.1, 0.25, 0.5, 1.0, 2.0, 4.0) returned the same 34 pre-computed Seurat clusters, making the cluster granularity sweep useless.  
**Cause**: The fallback candidate list includes `"seurat_clusters"` and any column matching `"cluster"` in the name. The pathogen Seurat object has such a column and it is found first, regardless of the requested `resolution`.  
**Fix**: Only reuse an existing column if its name explicitly encodes the requested resolution in the standard Seurat format (e.g. `RNA_snn_res.{resolution}` or `SCT_snn_res.{resolution}`). Remove `"seurat_clusters"` and the generic `"cluster"` pattern match from the fallback candidates entirely. When no resolution-specific column is found, always recompute via kNN + Louvain.

```r
# Current (broken): checks seurat_clusters and generic "cluster" pattern
candidates <- c(
  paste0("RNA_snn_res.", resolution),
  paste0("SCT_snn_res.", resolution),
  paste0("wsnn_res.", resolution),
  "seurat_clusters"              # ← remove this
)
# Also remove the grep("cluster", ...) fallback

# Fixed: only accept resolution-specific column names
candidates <- c(
  paste0("RNA_snn_res.", resolution),
  paste0("SCT_snn_res.", resolution),
  paste0("wsnn_res.", resolution)
)
# If none found → always recompute via kNN + Louvain
```

### Bug 2 — eval_heldout_predictivity O(n²) blowup (FIXED, commit 9c7684a)

**Was**: building a full n_genes × n_genes weight matrix (e.g. 11010 × 11010 = 970 MB) per CV fold, causing extreme GC pressure. For n_genes = 11010 × 5 folds, this took 2+ hours and produced astronomically wrong R² values (−10²⁹).  
**Fix applied**: when n_genes > 500, subsample 500 genes (seed = 42) for the R² estimate. Replaced W_abs matrix with a per-gene prediction loop + `rm(cor_train); gc()` between folds. Added `pmax(-1, pmin(1, r2))` clip. See commit `9c7684a`.  
**Residual note**: the 500-gene subsample is a heuristic. For datasets where n_genes < 500 (after filtering), the full gene set is used. The metric should be interpreted as an estimate of predictivity on a random gene subsample, not the full transcriptome.

---

## 5. Open FLAGs and status

| FLAG | Short description | Status |
|---|---|---|
| FLAG-01 | n_iter default: corrected to formula round(p(p-1)/39980) | **CORRECTED** (Phase 2c) |
| FLAG-02 | modulePreservation fallback to meancor z-score | **IMPLEMENTED** (schema captures preservation_method column) |
| FLAG-03 | Full R_score pair table must be saved (not filtered) | **IMPLEMENTED** (enforced in RobustnessResult schema; pair_scores_full.csv written) |
| FLAG-04 | Mixed gene IDs in casual GGM output | **IMPLEMENTED** (adapter maps symbols → AT-IDs at input) |
| FLAG-05 | Conditions pooled in casual GGM run — BLOCKER | **OPEN** (per-condition rerun planned; casual artifacts not usable as pipeline core) |
| FLAG-06 | HVG filter in casual GGM run — BLOCKER | **OPEN** (full-universe rerun on HPC required) |
| FLAG-07 | GOI lookup partners: long-format companion needed | **OPEN** (goi_partners_long.csv not yet implemented in goi_lookup.R) |
| FLAG-08 | Dev atlas Seurat objects archived off-machine | **OPEN** (not blocking; revisit when objects restored) |
| FLAG-09 | min_cells/coex_cutoff semantics conflated | **CORRECTED** (Phase 2c; coex_cutoff now filters on co-detection cell count) |
| FLAG-10 | Sign handling: positive-only vs \|pcor\| filter | **IMPLEMENTED** (keep_negative=FALSE default matches paper; toggle available) |
| FLAG-11 | WGCNA soft-thresholding on GGM networks | **RESOLVED** (dual-method: WGCNA p=1 + Louvain on large/small graphs) |
| FLAG-12 | R_score buries condition-specific pairs | **IMPLEMENTED** (characterize_condition_pattern() added; pair_condition_patterns.csv; 1.4M pairs) |
| FLAG-13 | GGM misses inducible/cell-specific TFs; pseudobulk captures them | **DOCUMENTED** (guidance added to ARCHITECTURE.md; both modes standard output) |
| FLAG-14 | Observation-point design as first-class variable | **STAGE 0-2 IMPLEMENTED** (generators, eval harness, normalization decided; cluster sweep bug open; metacell sweep incomplete) |

### FLAG-14 sub-items remaining open

- **Bug 1** (obs_cluster resolution): must fix before cluster granularity sweep is meaningful.
- **Metacell sweep**: 4 designs (t200/t100/t50/t25) not completed due to runtime. Needs PCA speedup (irlba) or reduced split-half reps.
- **Post-hoc BON3/WRKY sanity**: not run (sweep terminated before reaching post-hoc section). Run manually on the subcluster design.
- **Comparison methods**: hdWGCNA / CS-CORE / SuperCell / SEACells interface ready but not implemented.

---

## Session 2026-06-15

### Completed this session
- Bug #1 fixed: obs_cluster resolution fallback — removed seurat_clusters and
  grep("cluster") from candidate list; recompute is now default when no
  resolution-specific column exists; cluster_col populated correctly.
  Also fixed: cluster_col was NULL after recompute (shadowed loop variable).
  Regression test added. Committed and pushed.

- Phase 0 GGM review checklist completed (Axis 1–4, verdict: rerun required).
  Key blockers: 1.16 (pooled conditions), 1.9 (HVG-only), 2.8 (mixed gene IDs).
  Committed and pushed.

- Pseudobulk pipeline run (Steps 1–3 complete):
  zscore_gene + Spearman + obs_subcluster(298 pts)
  → pair_scores_full.csv (54M pairs), robustness_result.rds,
    pair_condition_patterns.csv
  All saved under results/pathogen_multiome/pseudobulk_zscore_spearman/

- Stage 3 (edge-threshold selection) established as new pipeline step,
  co-equal with Stage 1 (normalization) and Stage 2 (obs-point design).
  Same prior-free principle: stability-richness Pareto front.

- Stage 3 Phase 1 complete: 9 networks built and characterised.
  density_table.csv written. All 9 points valid.

### Key findings and decisions
- R_score cannot serve as an edge-density lever for pseudobulk Spearman
  networks. It is a discrete 4-value count of conditions (0.25/0.50/0.75/1.00),
  not a correlation-strength filter. Even R_score=1.00 leaves 42.5% density.
  In GGM networks it worked incidentally because pcor cutoff pre-sparsified the
  network; Spearman pseudobulk has no such pre-sparsification.
  → Correct lever: mean|r| = tanh(z_bar) (global threshold) OR per-gene top-k.

- Stage 3 Phase 2 (evaluation) ran but is INVALID. The stage3_eval_* adapters
  used 4 per-condition cor matrices as the observation axis:
    eff_rank saturated at ~2.8 (ceiling = 4 conditions)
    heldout_r2 degenerated (2-fold LOO)
    splithalf measured condition-structure divergence, not network stability
  → Must rewrite all adapters to use cell→298-obs-point axis (same as Stage 1/2).

- Evaluation adapter design decisions (confirmed):
    splithalf: cell 2-split → obs-point rebuild → threshold → Jaccard (primary)
               + Pearson of full upper triangle (secondary). 5 reps → 3 if >25min.
    eff_rank:  SVD of genes×298 obs-point matrix after threshold masking.
               (Was genes×4 fingerprint — wrong. Must be genes×298.)
    heldout:   obs-point 5-fold CV guilt-by-association R².
               (Was gene-hold-out LOO on 4 conditions — wrong.)
    null_gap:  real vs permuted edges above threshold. n_perm=10. Keep as-is.
    visible:   genes with ≥1 retained edge. Keep as-is.

- WGCNA is viable at appropriate density: power=1 gave 22 modules (grey 7.5%)
  at |r|≥0.42 in 3.6 min. Not viable at full density (74%: collapses to 1
  module). pickSoftThreshold deferred to after threshold selection.

- Environment trap (aarch64-darwin): data.table GForce segfaults.
  Fix: base-R + fread(nThread=1L). Multiple commits already in main.

### Open items (updated)
- FLAG-14 partially open:
    Stage 3 Phase 2: adapter rewrite + rerun NEXT (this handoff)
    metacell sweep: still pending (never reached in this session)
- Threshold confirmed as NOT determined yet (Phase 2 invalid)
-暫定モジュール results/pathogen_multiome/modules_absr042/ at |r|≥0.42:
  Louvain 6 modules, WGCNA-test 22 modules. Treat as preliminary only.
  Rebuild after threshold selection.
- pickSoftThreshold for WGCNA: deferred to after threshold confirmed.
- GGM rerun (per-condition, full gene universe): still pending (RERUN_PLAN.md).

### Stage 3 final decision (added post-sweep)
Edge threshold confirmed: global |r| = tanh(z_bar) ≥ 0.42

Prior-free justification (obs-point axis, 298 pseudobulk profiles):
  splithalf_jaccard = 0.9509
  eff_rank          = 116.87  (out of max 298)
  heldout_r2        = 0.5500
  null_gap          = 7.2M    (real_frac≈0.96, perm_frac≈0.00)
  visible_genes     = 5,450

Lever B (per-gene top-k) rejected: null_gap ≈ 1.19 across all k
(only 19% above random; eff_rank advantage was an artifact of
forced full-gene coverage, not network quality).

Within Lever A, |r|≥0.42 selected over stricter thresholds:
  - splithalf plateau above 0.42 (Δ=0.008 from 0.42→0.50)
  - visible_genes drops sharply above 0.46 (−26% at 0.50)
  - diminishing returns in stability, meaningful cost in coverage

Next steps:
  - Rebuild official modules at |r|≥0.42:
      WGCNA with pickSoftThreshold (deferred from session)
      Louvain (replace modules_absr042 with official run)
  - metacell sweep (FLAG-14, still pending)
  - GGM rerun per-condition full gene universe (RERUN_PLAN.md)
