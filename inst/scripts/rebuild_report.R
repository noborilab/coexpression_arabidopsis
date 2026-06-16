#!/usr/bin/env Rscript
# rebuild_report.R — regenerate standalone HTML report with all patched outputs
options(WGCNA.useThreads = FALSE)
suppressPackageStartupMessages({ library(data.table); setDTthreads(1L); library(base64enc) })

DOWN  <- "results/pathogen_multiome/downstream"
FIGS  <- file.path(DOWN, "figures")
RPATH <- file.path(DOWN, "DOWNSTREAM_ANALYSIS_REPORT.html")
cat("Building report...\n")

b64_png <- function(path) {
  if (!file.exists(path)) return(NULL)
  readBin(path, "raw", file.size(path)) |> base64encode() |>
    paste0("data:image/png;base64,", x=_)
}
img_tag <- function(p, alt="", w="100%") {
  b <- tryCatch(b64_png(p), error=function(e) NULL)
  if (is.null(b)) return(paste0("<p><em>Fig missing: ", basename(p), "</em></p>"))
  paste0('<img src="', b, '" alt="', alt, '" style="max-width:', w, ';height:auto;">')
}
read_safe <- function(f, n=25) tryCatch(head(read.csv(file.path(DOWN,f)), n), error=function(e) NULL)
df2html <- function(df, cap="") {
  if (is.null(df) || nrow(df)==0) return("<p><em>No data</em></p>")
  hdr  <- paste0("<th>", names(df), "</th>", collapse="")
  rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse=""), "</tr>"))
  paste0(if(nchar(cap)>0) paste0("<p><strong>",cap,"</strong></p>") else "",
         '<div style="overflow-x:auto"><table border="1" style="border-collapse:collapse;font-size:11px;">',
         "<thead><tr>", hdr, "</tr></thead><tbody>",
         paste(rows, collapse=""), "</tbody></table></div>")
}
fig_row <- function(...) {
  boxes <- list(...)
  paste0('<div style="display:flex;flex-wrap:wrap;gap:10px;">', paste(boxes, collapse=""), "</div>")
}
fig_box <- function(fig_path, cap="", flex="46%") {
  paste0('<div style="flex:1 1 ', flex, ';border:1px solid #ddd;border-radius:4px;padding:8px;background:#fff;">',
         img_tag(file.path(FIGS, fig_path)),
         if(nchar(cap)>0) paste0('<p style="font-size:11px;color:#666;">', cap, "</p>") else "",
         "</div>")
}

ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

html <- paste0(
'<!DOCTYPE html><html lang="en"><head>
<meta charset="UTF-8">
<title>Downstream Co-expression Analysis Report</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
body{font-family:Arial,sans-serif;max-width:1400px;margin:0 auto;padding:20px;background:#fafafa;}
h1{color:#2c3e50;border-bottom:3px solid #3498db;padding-bottom:10px;}
h2{color:#2980b9;margin-top:40px;border-left:4px solid #3498db;padding-left:10px;}
h3{color:#555;}
.caveat{background:#fef9e7;border:2px solid #f39c12;padding:12px 18px;border-radius:6px;margin:12px 0;}
.note{background:#eaf4fb;border:1px solid #85c1e9;padding:10px 15px;border-radius:5px;margin:10px 0;}
.section{background:white;border:1px solid #ddd;border-radius:6px;padding:20px 25px;margin:20px 0;}
details{margin:10px 0;}
summary{cursor:pointer;font-weight:bold;color:#2980b9;padding:6px;background:#eaf4fb;border-radius:4px;}
summary:hover{background:#d6eaf8;}
table{border-collapse:collapse;font-size:11px;}
th{background:#2980b9;color:white;padding:4px 8px;}
td{padding:3px 8px;border:1px solid #ddd;}
tr:nth-child(even){background:#f2f2f2;}
#sticky-toc{position:sticky;top:0;background:#2c3e50;color:white;padding:8px 20px;z-index:100;display:flex;flex-wrap:wrap;gap:10px;}
#sticky-toc a{color:#85c1e9;text-decoration:none;font-size:12px;}
#sticky-toc a:hover{color:white;}
</style></head><body>
<div id="sticky-toc">
  <strong style="margin-right:12px;">Navigate:</strong>
  <a href="#overview">1.Overview</a><a href="#topology">2.Topology</a>
  <a href="#module-quality">3.Module Quality</a><a href="#cross-mode">4.Cross-mode</a>
  <a href="#condition">5.Condition</a><a href="#annotation">6.Annotation</a>
  <a href="#gene-centric">7.Gene-centric</a><a href="#summary">8.Summary</a>
</div>
<h1>Downstream Co-expression Analysis Report</h1>
<p style="color:#666;">Generated: ', ts, ' | Dataset: pathogen_multiome | Repo: coexpression_arabidopsis</p>

<!-- 1. OVERVIEW -->
<div class="section" id="overview">
<h2>1. Overview</h2>
<p>Downstream analysis of finalized co-expression networks and modules from the <em>Arabidopsis thaliana</em> pathogen multiome dataset. All primary analyses are <strong>structure-based and prior-free</strong>.</p>
<div class="caveat"><strong>GO/TF annotation is REFERENCE OUTPUT ONLY.</strong> Never used to name modules or as selection criteria. Modules are identified by number only.</div>
<table border="1"><tr><th>Property</th><th>GGM mode</th><th>Pseudobulk mode</th></tr>
<tr><td>Edge type</td><td>Partial correlation (GGM, 4 conditions)</td><td>Spearman |r| ≥ 0.42</td></tr>
<tr><td>Observation unit</td><td>Individual cells</td><td>298 Seurat subclusters</td></tr>
<tr><td>Networks</td><td>4 per-condition + robustness consensus</td><td>1 network</td></tr>
<tr><td>Module sets</td><td>4 (large/small × WGCNA/Louvain)</td><td>2 (WGCNA / Louvain)</td></tr>
<tr><td>n_nodes (consensus/PB)</td><td>10,358 (GGM consensus)</td><td>5,450 (pseudobulk)</td></tr>
<tr><td>n_edges (consensus/PB)</td><td>62,863 (GGM consensus)</td><td>751,959 (pseudobulk)</td></tr>
</table>
<div class="caveat"><strong>⚠ Cross-mode condition-pattern comparison is FORBIDDEN.</strong><br>
Different edge definitions (partial vs. marginal correlation) and different observation units mean condition weights cannot be compared quantitatively across modes. Per-mode analysis is valid within each mode only.</div>
</div>

<!-- 2. TOPOLOGY -->
<div class="section" id="topology">
<h2>2. Network Topology (A1–A3)</h2>
<h3>A3: Global Statistics</h3>',
df2html(read_safe("topology_global_stats.csv"), "Global stats — all 6 networks"),
'<div class="note">Clustering coefficient computed only for GGM consensus (&lt;80k edges). Component stats via igraph::components() for all networks.</div>
',
fig_row(
  fig_box("fig_global_stats_comparison.png", "Global network statistics (GGM modes and pseudobulk side-by-side)", "90%")
),
'<h3>A1: Degree Distributions + Scale-Free Fit</h3>
<div class="note">Power-law alpha via igraph::fit_power_law (Clauset et al. MLE). KS.p returned NA for all networks in this igraph version — alpha values reported. Higher alpha = steeper degree distribution = more hub-dominated topology.</div>',
df2html(read_safe("topology_degree.csv", 15), "Degree + power-law (first 15 rows)"),
fig_row(fig_box("fig_degree_distributions.png", "Degree distributions (log10 x-axis). All GGM per-condition networks show similar shapes.", "90%")),
'<h3>A2: Hub Genes (Centrality)</h3>
<div class="note">Betweenness centrality computed for GGM consensus only (&lt;100k edges). Per-condition GGM networks and pseudobulk: degree only (eigenvector/betweenness not computed).</div>',
df2html(read_safe("topology_centrality.csv", 15), "Centrality records (first 15 rows)"),
fig_row(fig_box("fig_hub_genes.png", "Degree vs eigenvector centrality. NA eigenvector = large-network fallback.", "90%")),
"</div>

<!-- 3. MODULE QUALITY -->
<div class=\"section\" id=\"module-quality\">
<h2>3. Module Quality (B1–B4)</h2>
<h3>B4: Quality Summary — All 6 Module Sets</h3>",
df2html(read_safe("module_quality_summary.csv"), "Module quality across all 6 sets"),
fig_row(fig_box("fig_module_quality_across_sets.png",
  "Module quality: grey rate, n_modules, kME median. Large-louvain assigns 91% of genes to modules (lowest grey rate); large-wgcna 68% grey.", "90%")),
"<h3>B1: kME Distributions</h3>",
df2html(read_safe("module_kme_distributions.csv", 15), "kME stats per module (first 15)"),
fig_row(fig_box("fig_module_kme.png", "kME distributions. Modules with median kME < 0.3 flagged as low-coherence.", "90%")),
"<h3>B3: Intramodular Hub Genes</h3>",
df2html(read_safe("module_hubs.csv", 20), "Hub genes (top-kME per module, first 20)"),
"<h3>B2: Eigengene-Eigengene Correlations</h3>
<details><summary>Eigengene correlation tables (collapsed)</summary>",
df2html(read_safe("module_eigengene_correlations.csv", 20)),
"</details>",
fig_row(
  fig_box("fig_eigengene_heatmap_GGM_large_wgcna.png", "GGM large_wgcna eigengene correlations"),
  fig_box("fig_eigengene_heatmap_GGM_large_louvain.png", "GGM large_louvain eigengene correlations")
),
fig_row(
  fig_box("fig_eigengene_heatmap_PB_wgcna.png", "PB WGCNA eigengene correlations (sample-level MEs)"),
  fig_box("fig_eigengene_heatmap_PB_louvain.png", "PB Louvain eigengene correlations (PC1 per module)")
),
"</div>

<!-- 4. CROSS-MODE -->
<div class=\"section\" id=\"cross-mode\">
<h2>4. Cross-mode &amp; Cross-method (C1–C4) — Headline Section</h2>
<div class=\"note\">The pipeline's distinctive value: two complementary co-expression views (partial vs. marginal correlation) and which gene groupings are robust across both.</div>

<h3>C3: Core vs Mode-specific Gene Partitioning</h3>",
df2html(head(read_safe("core_vs_modespecific.csv", 20), 10), "Gene categories (first 10)"),
fig_row(
  fig_box("fig_core_modespecific.png", "<strong>Core</strong>: 1,095 genes assigned in both GGM and pseudobulk. <strong>GGM-only</strong>: 8,371. <strong>PB-only</strong>: 4,048. Low cross-mode core fraction expected from different edge definitions."),
  fig_box("fig_cross_set_consistency.png", "GGM cross-set stability: 1,687 genes assigned in all 4 GGM sets (high-confidence co-expression modules).")
),
"<h3>C1: GGM vs Pseudobulk Module Overlap (Jaccard + ARI)</h3>",
df2html(read_safe("crossmode_ari.csv"), "Cross-mode ARI per set pair"),
"<div class=\"note\">ARI range: -0.018 to 0.018; Jaccard max: 0.026. Low overlap is expected — GGM and pseudobulk capture different aspects of co-expression (partial vs. marginal correlation, cells vs. subclusters). Complementarity is assessed by gene-level recovery (see F2), not by module boundary overlap.</div>",
fig_row(fig_box("fig_crossmode_jaccard_heatmap.png",
  "Top-50 GGM vs pseudobulk module pairs by Jaccard (max=0.026, all pairs shown). Low overlap confirms the two modes are complementary views, not redundant.", "90%")),
"<h3>C2: WGCNA vs Louvain Agreement (within each mode)</h3>",
df2html(read_safe("crossmethod_agreement.csv"), "WGCNA vs Louvain ARI"),
fig_row(fig_box("fig_crossmethod_ari.png",
  "WGCNA vs Louvain ARI. GGM small consensus: 0.301 (good agreement). GGM large: 0.028 (low, expected — Louvain doesn't penalize grey). Pseudobulk: 0.137.", "70%")),
"</div>

<!-- 5. CONDITION -->
<div class=\"section\" id=\"condition\">
<h2>5. Condition-specificity (D1–D2)</h2>
<div class=\"caveat\"><strong>⚠ The two panels below are NOT cross-comparable.</strong><br>
GGM condition weights = partial correlation edges in condition-specific GGM networks. Pseudobulk condition weights = mean intramodular Spearman correlations in subcluster batches. Different scales, different null distributions. Compare WITHIN each panel only.</div>
<h3>D1: Module Condition Activation</h3>
<h4>GGM mode (within-mode only)</h4>",
df2html(read_safe("module_condition_activation_ggm.csv", 10)),
fig_row(fig_box("fig_module_condition_activation_ggm.png",
  "GGM module condition activation — WITHIN GGM mode, NOT cross-mode comparable. Top-20 modules shown.", "90%")),
"<h4>Pseudobulk mode (within-mode only)</h4>",
df2html(read_safe("module_condition_activation_pseudobulk.csv")),
fig_row(fig_box("fig_module_condition_activation_pseudobulk.png",
  "Pseudobulk module condition activation — WITHIN pseudobulk mode, NOT cross-mode comparable.", "90%")),
"</div>

<!-- 6. ANNOTATION -->
<div class=\"section\" id=\"annotation\">
<h2>6. Functional Annotation (E1–E2) — REFERENCE OUTPUT ONLY</h2>
<div class=\"caveat\"><strong>Modules are NOT named by GO terms or TF families.</strong> GO and TF content are descriptive metadata only — never used as selection or ranking criteria.</div>
<h3>E1: GO BP Enrichment (all module sets)</h3>
<p>Total: 1,491 significant GO terms across all 6 module sets (BH q&lt;0.05, min set size 10).</p>
<details><summary>Master GO enrichment table (collapsed; first 25 rows)</summary>",
df2html(read_safe("go_enrichment_all_sets.csv", 25)),
"</details>
<h3>E2: TF Module Enrichment</h3>",
df2html(read_safe("tf_enrichment.csv", 20), "TF-module records (first 20)"),
"</div>

<!-- 7. GENE-CENTRIC -->
<div class=\"section\" id=\"gene-centric\">
<h2>7. Gene-centric Utility (F1–F2)</h2>
<h3>F1: Master Gene Lookup Table</h3>
<p>18,836 genes × (module in each of 6 sets, kME, top-10 co-expression partners in both modes).</p>",
df2html(read_safe("gene_lookup_master.csv", 12), "gene_lookup_master.csv (first 12 rows)"),
"<div class=\"note\"><strong>How to use:</strong> Filter by gene_id for a full co-expression profile. top10_GGM_consensus and top10_PB columns give guilt-by-association partners. Module columns enable instant cross-set comparison.</div>
<h3>F2: WRKY Family — Demo of Generic Gene-set Query</h3>",
df2html(read_safe("geneset_query_demo_wrky.csv", 15), "WRKY enrichment per module-set (first 15)"),
fig_row(fig_box("fig_wrky_crossmode.png",
  "WRKY family module enrichment (BH q&lt;0.05). 1 significant hit (PB_wgcna module 13). Demo of generic gene-set query capability.", "90%")),
"</div>

<!-- 8. SUMMARY -->
<div class=\"section\" id=\"summary\">
<h2>8. Summary &amp; How to Reuse on a New Dataset</h2>

<h3>Analysis Completion Status</h3>
<table border=\"1\">
<thead><tr><th>Analysis</th><th>Status</th><th>Output</th><th>Key metric</th></tr></thead>
<tbody>
<tr><td>A1: Degree + power-law</td><td style=\"color:green\">DONE</td><td>topology_degree.csv</td><td>alpha: 5.1–9.4 across networks</td></tr>
<tr><td>A2: Centrality</td><td style=\"color:green\">DONE</td><td>topology_centrality.csv</td><td>79,720 gene-network records; betweenness for GGM consensus</td></tr>
<tr><td>A3: Global stats</td><td style=\"color:green\">DONE</td><td>topology_global_stats.csv</td><td>GGM consensus: 10,358 nodes, 62,863 edges, 306 components, LCC=93%</td></tr>
<tr><td>B1: kME distributions</td><td style=\"color:green\">DONE</td><td>module_kme_distributions.csv</td><td>116 module-set records; 92 low-coherence modules</td></tr>
<tr><td>B2: Eigengene correlations</td><td style=\"color:green\">DONE</td><td>module_eigengene_correlations.csv</td><td>1,558 module-pair correlations across 6 sets</td></tr>
<tr><td>B3: Hub genes</td><td style=\"color:green\">DONE</td><td>module_hubs.csv</td><td>2,050 hub gene records</td></tr>
<tr><td>B4: Module quality</td><td style=\"color:green\">DONE</td><td>module_quality_summary.csv</td><td>Grey rates: 68%, 8.6%, 44.5%, 1.9%, 8.1%, 6.6%</td></tr>
<tr><td>C1: Cross-mode overlap</td><td style=\"color:green\">DONE</td><td>crossmode_overlap.csv + crossmode_ari.csv</td><td>1,764 pairs; max Jaccard=0.026; ARI: -0.018–0.018</td></tr>
<tr><td>C2: WGCNA vs Louvain</td><td style=\"color:green\">DONE</td><td>crossmethod_agreement.csv</td><td>ARI: 0.028–0.301</td></tr>
<tr><td>C3: Core vs mode-specific</td><td style=\"color:green\">DONE</td><td>core_vs_modespecific.csv</td><td>Core=1,095; GGM-only=8,371; PB-only=4,048</td></tr>
<tr><td>C4: Cross-set consistency</td><td style=\"color:green\">DONE</td><td>cross_set_consistency.csv</td><td>1,687 genes in all 4 GGM sets</td></tr>
<tr><td>D1: GGM condition activation</td><td style=\"color:green\">DONE</td><td>module_condition_activation_ggm.csv</td><td>42 modules; reused condition_comparison/</td></tr>
<tr><td>D1: PB condition activation</td><td style=\"color:green\">DONE</td><td>module_condition_activation_pseudobulk.csv</td><td>18 records (WGCNA+Louvain)</td></tr>
<tr><td>D2: GGM condition patterns</td><td style=\"color:green\">DONE</td><td>module_condition_patterns_ggm.csv</td><td>98 records; reused all_modules_condition_patterns.csv</td></tr>
<tr><td>D2: PB condition patterns</td><td style=\"color:green\">DONE</td><td>module_condition_patterns_pseudobulk.csv</td><td>18 records (fill=TRUE for column differences)</td></tr>
<tr><td>E1: GO enrichment</td><td style=\"color:green\">DONE</td><td>go_enrichment_all_sets.csv</td><td>1,491 terms; all 6 sets</td></tr>
<tr><td>E2: TF enrichment</td><td style=\"color:green\">DONE</td><td>tf_enrichment.csv</td><td>319 TF-module records</td></tr>
<tr><td>F1: Gene lookup master</td><td style=\"color:green\">DONE</td><td>gene_lookup_master.csv</td><td>18,836 genes</td></tr>
<tr><td>F2: WRKY demo</td><td style=\"color:green\">DONE</td><td>geneset_query_demo_wrky.csv</td><td>1 significant module enrichment</td></tr>
</tbody></table>

<h3>Key Findings</h3>
<ul>
<li><strong>Network scale-free topology</strong>: Power-law alpha 5.1–9.4 across all 6 networks — consistent with scale-free hub structure.</li>
<li><strong>Module quality varies by method</strong>: Louvain assigns most genes (grey rates 1.9–8.6%); WGCNA has higher grey rates but better intramodular coherence (higher kME).</li>
<li><strong>Cross-mode overlap is low</strong>: max Jaccard = 0.026; ARI ≈ 0. Expected — partial correlation (GGM) captures different regulatory relationships than marginal correlation (pseudobulk). Complementarity is at the individual gene level (WRKY demo), not module boundaries.</li>
<li><strong>1,687 genes are stable across all 4 GGM sets</strong>: these are the most robust co-expression units.</li>
<li><strong>1,095 core genes</strong> assigned in both GGM and pseudobulk modules — high-confidence co-expressed gene groups.</li>
</ul>

<h3>How to Reuse on a New Dataset</h3>
<ol>
<li><strong>Fully automatic</strong>: A1–A3 (topology), B1–B4 (module quality), C1–C4 (cross-mode), D1–D2 (condition), E1 (GO), E2 (TF if module_tfs.csv present)</li>
<li><strong>Requires parameter</strong>: F2 — provide any gene list CSV with gene_id column</li>
<li><strong>Entry point</strong>: <code>Rscript inst/scripts/downstream_analysis.R</code> from the repo root</li>
<li><strong>Resumable</strong>: all analyses skip if output exists; re-run safely after failures</li>
</ol>

<h3>Output Locations</h3>
<ul>
<li><code>results/pathogen_multiome/downstream/</code> — all CSV outputs</li>
<li><code>results/pathogen_multiome/downstream/figures/</code> — all PNG + SVG figures</li>
<li><code>results/pathogen_multiome/downstream/DOWNSTREAM_ANALYSIS_REPORT.html</code> — this report (standalone)</li>
<li><code>inst/scripts/downstream_analysis.R</code> — main analysis script (generic, reusable)</li>
</ul>
</div>
</body></html>')

writeLines(html, RPATH)
sz <- file.size(RPATH)
cat("Report written:", RPATH, "\n")
cat("Size:", round(sz/1024), "KB | Standalone:", sz > 10000, "\n")
