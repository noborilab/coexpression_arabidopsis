# Extended Gene Coexpression Analysis Pipeline — Background & Design Substrate

This document is the context handoff for a new chat dedicated to designing and building a
reusable gene co-expression analysis pipeline (working name: **Extended Gene Coexpression
Analysis Pipeline**) as a lab tool / GitHub repo (`noborilab/<name>`). It captures what was
learned in a prior project so the new chat starts fully informed. It is a strategy substrate,
not a locked spec.

-----

## 1. Goal of the pipeline

Turn the ad-hoc co-expression analyses already done (hard-coded to two datasets) into a
**reusable pipeline that takes any 10x Genomics Seurat object and produces tissue/context-robust
co-expression modules with biological interpretation**. Scope is explicitly **beyond the CZL
collaboration** — it should apply to any of the lab’s 10x datasets (PRIMER/bystander immunity
data, future atlases, csRNA-seq-linked work, etc.).

Two analysis modes must BOTH be first-class citizens:

- **Pseudobulk mode** — cluster/sample pseudobulk + marginal correlation (what was done on the
  dev atlas and pathogen data).
- **SingleCellGGM mode** — cell-level graphical Gaussian model / partial correlation (already
  run casually on the pathogen data; needs review).

**cross-stratum robustness** (reproducibility of co-expression across independent strata —
organs, conditions, samples) is a valuable layer where the data supports it (depends on the
input dataset’s structure), and should be an optional/configurable layer, not hard-wired.

-----

## 2. Where this came from (the prior project, condensed)

Originally a ligand–receptor (LR) discovery effort for the CZL collaboration (Zipfel lab),
analysing propeptide/receptor co-expression on two Arabidopsis single-nucleus datasets:

- **Dev atlas**: pseudobulk count tables, ~34k genes, clusters nested in 10 organs/stages
  (seed_3d/425d, seedling_6d/9d/15d, flower, rosette_21d/30d, stem, silique). Pseudobulk only;
  the cell-level Seurat objects are archived off-machine (large) → SingleCellGGM on dev is
  deferred until they’re restored.
- **Pathogen multiome**: snRNA+ATAC of Arabidopsis leaf, conditions Mock/DC3000/AvrRpt2/AvrRpm1,
  65k nuclei, 15 samples, ~298–428 sub-clusters. RNA pseudobulk precomputed. SingleCellGGM was
  run on this (cell-level) casually by Tatsuya — **needs review before being used as a pipeline
  core**.

**Key scientific finding that reframed everything:** pseudobulk co-expression cannot recover
classical LR pairs because plant peptide–receptor signalling is overwhelmingly **paracrine**
(ligand and receptor in different cell types → zero/negative co-expression). This is a hard
limit of ALL co-expression methods (single-cell GGM included; single-cell resolution makes
paracrine pairs HARDER, not easier). So the project pivoted away from “LR pair discovery” to
**co-expression module discovery + context interpretation** (grouping genes by expression
context and interpreting that context), delivered as input/hypotheses for CZL experiments — not
a standalone computational paper.

-----

## 3. Methods already built (reusable parts to migrate into the pipeline)

All in R, pseudobulk-based, under a prior `out/CZL_discovery/` tree. The reusable machinery:

- **Robustness statistic (`robuststat` library).** Per-stratum Spearman via
  rank-transform-then-Pearson; **fixed-evidence indicator** `I_s = 1[z_s ≥ k·SE_s]` (NOT a fixed
  ρ cutoff — a small stratum needs a larger ρ; `z_s=atanh(ρ_s)`, `SE_s=1/sqrt(n_s−3)`, `k≈1.64`,
  calibrated on positive controls); aggregated `R_score = Σ w_s·I_s / Σ w_s`,
  `w_s=sqrt(min(n_s,30)−3)`. Null: **analytic weighted Poisson-binomial** (estimate per-stratum
  null indicator prob π_s by pooling matched-permutation draws across pairs in the same
  (n_s, expression×detection bin), then Poisson-binomial upper tail → resolved p, no Monte-Carlo
  floor) → BH-FDR. Companions: random-effects z̄, between-stratum τ², anchored conditional
  Spearman (for narrow-ligand asymmetry), co-detection Jaccard, promiscuity flag.
- **Cross-dataset / cross-context replication** as a bonus annotation (no AND-gate): a pair/edge
  robust in dev (organ-controlled) AND pathogen (condition-controlled) gets a “star”.
- **WGCNA module construction** (genome-wide). Signed network; soft-power by scale-free fit;
  module merge threshold tunes granularity. A granularity sweep over (power × merge) was used to
  pick the set.
- **Per-module context interpretation**: dev organ activity (eigengene), pathogen treatment
  activity + Δ(infected−Mock), GO BP enrichment, curated-set anchors, plain-language label.
- **Module preservation** across contexts (WGCNA modulePreservation Zsummary; with a lightweight
  mean-intramodular-|cor| z-score fallback when modulePreservation times out).
- **Hierarchy / drill-down**: top-level modules (coarse merge) with nested sub-modules (fine
  merge) from the same network, so large modules can be dissected.
- **Hub genes** (kME), **TF intersection** (regulator hints per module), **gene-of-interest
  lookup table** (per gene: module, kME, hub flag, preservation).

Important method-benchmark result (sets expectations): on GO co-functional pair recovery, the
tissue-robust R_score did **NOT** beat naive Spearman or WGCNA (AUPRC ~0.20 for all three). The
de-confounding removes organ-driven pairs that are modestly less co-functional (19.4% vs 24.3%),
but does not improve global recovery. **WGCNA beat hclust-on-organ-means on curated-set
fold-enrichment (50× vs 27×).** Takeaway for the pipeline: the value is producing clean,
interpretable, context-annotated modules — not a claim that robustness “beats” standard methods.

-----

## 4. The recommended pipeline architecture (discussed, not yet built)

The cleanest synthesis of the two modes:

> **per-context network estimation → cross-context preservation/robustness → modules +
> interpretation**

- **Input**: a Seurat object (10x) + a specification of the stratum/context variable(s)
  (organ, condition, sample, timepoint…) and the cell-grouping for pseudobulk.
- **Estimation mode (configurable):**
  - *Pseudobulk*: aggregate to cluster/sample pseudobulk, marginal Spearman (+ robustness
    statistic per stratum).
  - *SingleCellGGM*: cell-level partial correlation (graphical Gaussian model). Partial
    correlation removes indirect edges and partially absorbs the tissue-identity confound by
    conditioning on other genes — a genuine advantage over marginal correlation for network
    quality. Same fundamental paracrine limitation for LR pairs remains.
- **Robustness layer (optional, data-dependent):** cross-stratum reproducibility (R_score) and/or
  cross-dataset replication — only meaningful when the input has multiple independent strata.
- **Output layer (shared by both modes):** WGCNA-style modules (or GGM modules), preservation,
  hierarchy/hubs, GO + TF + curated-anchor interpretation, gene-of-interest lookup.
- **Domain plugins:** keep generic core; push dataset-specific pieces (e.g. CZL gene lists,
  custom curated sets) into swappable config/plugins.

Practical constraints/notes:

- SingleCellGGM needs cell-level matrices and is the heavier compute; pseudobulk is light.
- WGCNA `modulePreservation` is slow on large cell-state matrices (timed out at 40 min on
  ~1500 samples × 11k genes) → cap/representative-subsample or use the fallback proxy.
- TF list that worked: lab’s own motif metadata
  `…/from_Ben/for_tatsuya/data/motifs-2026/Athaliana_motifs_metadata.tsv` (673 TFs, col1 =
  AT-ID `motif_id`, plus symbol/family/class) — better than PlantTFDB auto-download which failed.

-----

## 5. Suggested build sequence

- **Phase 0 — review the existing pathogen SingleCellGGM run** (done casually by Tatsuya).
  Check parameters (pcor cutoff, min-cells, subsampling iterations), output format, reproducibility,
  and whether it can serve as the pipeline core. This gates the core design.
- **Phase 1 — lock the core design**: input contract (Seurat object + stratum spec), estimation
  modes, optional robustness layer, shared output/interpretation layer, plugin boundary.
- **Phase 2 — repo (`noborilab/<name>`)**: CLI, config-driven, tests; migrate the robuststat +
  WGCNA + interpretation parts; both pseudobulk and SingleCellGGM modes wired to one output schema.

Reference method paper: SingleCellGGM — Xu, Wang & Ma (2024), *Cell Reports Methods* 4, 100813
(“SingleCellGGM enables gene expression program identification…”). Single-cell graphical Gaussian
model: iterative random 2,000-gene subsampling, takes the minimum |pcor| across iterations as the
conservative final pcor; retains pcor ≥ 0.03 in ≥10 cells; outputs gene expression programs (GEPs)
= co-expression modules. Cells are the observations (not metacells/pseudobulk).

-----

## 6. Lab tooling conventions (for the repo)

The lab already ships standalone tools as repos (SnailFinder, CRISPR guide designer,
nobori-lab-db, receipt-processor, meeting-log). Follow the same pattern: `noborilab/<name>`,
config-driven CLI, tests, README. Work is typically driven via Claude Code on the lab machine.

-----

## 7. What this new chat should produce

1. A review checklist / plan for the existing pathogen SingleCellGGM run (Phase 0).
1. A locked core architecture (Phase 1): input contract, the two estimation modes, the optional
   robustness layer, shared output schema, plugin boundary.
1. A repo scaffold + an executable Claude Code build spec (Phase 2).

Note: the CZL-specific deliverable (module results for the Zipfel lab) is being finished in a
separate chat and is essentially complete — this pipeline effort is the generalization, not the
CZL delivery.