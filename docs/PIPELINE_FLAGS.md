# Pipeline Design Flags & Decisions

Running log of issues found during Phase 0–1 that affect implementation or
interpretation. Update each flag's Status line as work progresses.

---

## FLAG-01: n_iter default
**Phase**: 1
**Issue**: The casual pathogen GGM run used ~20,000 iterations (formula-driven).
Pipeline spec sets n_iter = 100 (SingleCellGGM paper default).
**Decision**: Keep 100. The min-pcor aggregation converges well below 20k.
Verify convergence when Mock condition runs first in the rerun.
**Status**: Decided — no action needed in Phase 2.

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
