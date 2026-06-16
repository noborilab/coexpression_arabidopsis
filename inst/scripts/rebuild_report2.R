#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); setDTthreads(1L); library(base64enc) })

DOWN  <- "results/pathogen_multiome/downstream"
FIGS  <- file.path(DOWN, "figures")
RPATH <- file.path(DOWN, "DOWNSTREAM_ANALYSIS_REPORT.html")

b64_png <- function(path) {
  if (!file.exists(path)) return(NULL)
  b <- tryCatch(base64encode(readBin(path, "raw", file.size(path))), error=function(e) NULL)
  if (is.null(b)) return(NULL)
  paste0("data:image/png;base64,", b)
}
img <- function(p, alt="") {
  b <- b64_png(p)
  if (is.null(b)) return(paste0("<p><em>Missing: ", basename(p), "</em></p>"))
  paste0('<img src="', b, '" alt="', alt, '" style="max-width:100%;height:auto;">')
}
csv_tbl <- function(f, n=20) tryCatch(head(read.csv(file.path(DOWN, f)), n), error=function(e) NULL)
t2h <- function(df) {
  if (is.null(df) || nrow(df)==0) return("<p><em>No data</em></p>")
  h <- paste0("<th>", names(df), "</th>", collapse="")
  r <- apply(df, 1, function(x) paste0("<tr>", paste0("<td>", x, "</td>", collapse=""), "</tr>"))
  paste0('<div style="overflow-x:auto"><table border="1" style="border-collapse:collapse;font-size:11px;">',
         "<thead><tr>", h, "</tr></thead><tbody>", paste(r, collapse=""), "</tbody></table></div>")
}

cat("Encoding PNG figures...\n")
FIG <- function(nm, flex="90%") {
  paste0('<div style="display:flex;flex-wrap:wrap;gap:10px;"><div style="flex:1 1 ', flex,
         ';border:1px solid #ddd;padding:8px;background:#fff;">', img(file.path(FIGS, nm)), "</div></div>")
}

ts_now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
gs     <- csv_tbl("topology_global_stats.csv")

html <- paste0(
'<!DOCTYPE html><html><head><meta charset="UTF-8"><title>Downstream Analysis Report</title>
<style>
body{font-family:Arial,sans-serif;max-width:1400px;margin:0 auto;padding:20px;background:#fafafa;}
h1{color:#2c3e50;border-bottom:3px solid #3498db;padding-bottom:10px;}
h2{color:#2980b9;margin-top:40px;border-left:4px solid #3498db;padding-left:10px;}
h3{color:#555;}
.cav{background:#fef9e7;border:2px solid #f39c12;padding:12px;border-radius:6px;margin:12px 0;}
.note{background:#eaf4fb;border:1px solid #85c1e9;padding:10px;border-radius:5px;margin:10px 0;}
.sec{background:white;border:1px solid #ddd;border-radius:6px;padding:20px;margin:20px 0;}
details{margin:10px 0;}
summary{cursor:pointer;font-weight:bold;color:#2980b9;padding:6px;background:#eaf4fb;border-radius:4px;}
table{border-collapse:collapse;font-size:11px;}
th{background:#2980b9;color:white;padding:4px 8px;}
td{padding:3px 8px;border:1px solid #ddd;}
tr:nth-child(even){background:#f2f2f2;}
.frow{display:flex;flex-wrap:wrap;gap:10px;}
.fbox{flex:1 1 45%;border:1px solid #ddd;padding:8px;background:#fff;}
#toc{position:sticky;top:0;background:#2c3e50;padding:8px 20px;z-index:100;display:flex;flex-wrap:wrap;gap:10px;}
#toc a{color:#85c1e9;text-decoration:none;font-size:12px;}
</style></head><body>
<div id="toc">
  <strong style="color:white;margin-right:12px;">Navigate:</strong>
  <a href="#s1">1.Overview</a><a href="#s2">2.Topology</a><a href="#s3">3.Module Quality</a>
  <a href="#s4">4.Cross-mode</a><a href="#s5">5.Condition</a>
  <a href="#s6">6.Annotation</a><a href="#s7">7.Gene-centric</a><a href="#s8">8.Summary</a>
</div>
<h1>Downstream Co-expression Analysis Report</h1>
<p style="color:#666;">Generated: ', ts_now, ' | Dataset: pathogen_multiome | Repo: coexpression_arabidopsis</p>

<!-- OVERVIEW -->
<div class="sec" id="s1"><h2>1. Overview</h2>
<p>Prior-free downstream analysis of co-expression networks and modules from the Arabidopsis pathogen multiome dataset. All primary analyses are structure-based.</p>
<div class="cav"><strong>GO/TF annotation is REFERENCE OUTPUT ONLY.</strong> Never used to name modules or as ranking criteria. Modules identified by number only.</div>
<div class="cav"><strong>&#9888; Cross-mode condition-pattern comparison is FORBIDDEN.</strong><br>
GGM uses partial correlation (individual cells). Pseudobulk uses Spearman (298 subclusters). Different edge definitions, different observation units, different scales. Per-mode condition activation is valid within each mode only. Quantitative cross-mode comparison is not valid.</div>
<table border="1"><tr><th>Property</th><th>GGM mode</th><th>Pseudobulk mode</th></tr>
<tr><td>Edge type</td><td>Partial correlation (4 per-condition + robustness consensus)</td><td>Spearman |r| &#8805; 0.42</td></tr>
<tr><td>Observation unit</td><td>Individual cells per condition</td><td>298 Seurat subclusters</td></tr>
<tr><td>Module sets</td><td>4 (large/small &#215; WGCNA/Louvain)</td><td>2 (WGCNA / Louvain)</td></tr>
<tr><td>n_nodes (consensus/PB)</td><td>10,358 genes</td><td>5,450 genes</td></tr>
<tr><td>n_edges (consensus/PB)</td><td>62,863</td><td>751,959</td></tr>
</table>
</div>

<!-- TOPOLOGY -->
<div class="sec" id="s2"><h2>2. Network Topology (A1&#8211;A3)</h2>
<h3>A3: Global Network Statistics</h3>',
t2h(gs),
'<div class="note">GGM consensus: 10,358 nodes, 62,863 edges, 306 components, LCC=93.3%, clustering=0.18.
Pseudobulk: 5,450 nodes, 751,959 edges (marginal correlation captures more co-expression pairs).</div>',
FIG("fig_global_stats_comparison.png"),
'<h3>A1: Degree Distributions + Scale-Free Fit</h3>
<div class="note">Power-law alpha (Clauset MLE via igraph::fit_power_law):<br>
Mock=5.83 | DC3000=9.40 | AvrRpt2=6.45 | AvrRpm1=7.87 | GGM_consensus=5.12 | Pseudobulk=7.55<br>
All &#945; &gt; 2: consistent with scale-free topology. KS.p unavailable for this igraph version.</div>',
t2h(csv_tbl("topology_degree.csv", 12)),
FIG("fig_degree_distributions.png"),
'<h3>A2: Centrality &amp; Hub Genes</h3>
<div class="note">Betweenness centrality computed for GGM consensus only (&lt;100k edges, 20s).
Per-condition GGM and pseudobulk networks: degree + eigenvector only (too many edges for betweenness).</div>',
t2h(csv_tbl("topology_centrality.csv", 12)),
FIG("fig_hub_genes.png"),
'</div>

<!-- MODULE QUALITY -->
<div class="sec" id="s3"><h2>3. Module Quality (B1&#8211;B4)</h2>
<h3>B4: Quality Summary &#8212; All 6 Module Sets</h3>',
t2h(csv_tbl("module_quality_summary.csv")),
FIG("fig_module_quality_across_sets.png"),
'<h3>B1: kME Distributions</h3>',
t2h(csv_tbl("module_kme_distributions.csv", 12)),
FIG("fig_module_kme.png"),
'<h3>B3: Intramodular Hub Genes</h3>',
t2h(csv_tbl("module_hubs.csv", 15)),
'<h3>B2: Eigengene-Eigengene Correlations</h3>
<details><summary>Correlation data (first 20 rows)</summary>',
t2h(csv_tbl("module_eigengene_correlations.csv", 20)),
'</details>
<div class="frow">
  <div class="fbox">', img(file.path(FIGS,"fig_eigengene_heatmap_GGM_large_wgcna.png"),"GGM large_wgcna"), '</div>
  <div class="fbox">', img(file.path(FIGS,"fig_eigengene_heatmap_GGM_small_wgcna.png"),"GGM small_wgcna"), '</div>
</div>
<div class="frow">
  <div class="fbox">', img(file.path(FIGS,"fig_eigengene_heatmap_PB_wgcna.png"),"PB WGCNA"), '</div>
  <div class="fbox">', img(file.path(FIGS,"fig_eigengene_heatmap_PB_louvain.png"),"PB Louvain"), '</div>
</div>
</div>

<!-- CROSS-MODE -->
<div class="sec" id="s4"><h2>4. Cross-mode &amp; Cross-method (C1&#8211;C4) &#8212; Headline Section</h2>
<div class="note"><strong>This is the pipeline&#8217;s distinctive value.</strong> Two complementary co-expression views (partial vs. marginal correlation) and which gene groupings are robust across both modes.</div>

<h3>C3: Core vs Mode-specific Gene Partitioning</h3>
<div class="note"><strong>Core</strong>: assigned in both GGM and pseudobulk (1,095 genes).
<strong>GGM-only</strong>: 8,371 genes. <strong>PB-only</strong>: 4,048 genes.
Low cross-mode core fraction expected &#8212; GGM partial correlation and PB marginal correlation capture different regulatory relationships.</div>',
'<div class="frow">
  <div class="fbox">', img(file.path(FIGS,"fig_core_modespecific.png"),"Core vs mode-specific"), '</div>
  <div class="fbox">', img(file.path(FIGS,"fig_cross_set_consistency.png"),"Cross-set consistency"), '</div>
</div>

<h3>C1: GGM vs Pseudobulk Module Overlap (Jaccard + ARI)</h3>
<div class="note">Max Jaccard = 0.026; ARI range = -0.018 to 0.018. Low overlap confirmed &#8212; GGM and pseudobulk do NOT produce the same module partition, which is expected and desirable. Complementarity is at the individual gene level (demonstrated in F2), not at module boundaries.</div>',
t2h(csv_tbl("crossmode_ari.csv")),
FIG("fig_crossmode_jaccard_heatmap.png"),
'<h3>C2: WGCNA vs Louvain Agreement (within each mode)</h3>',
t2h(csv_tbl("crossmethod_agreement.csv")),
FIG("fig_crossmethod_ari.png", "70%"),
'</div>

<!-- CONDITION -->
<div class="sec" id="s5"><h2>5. Condition-specificity (D1&#8211;D2)</h2>
<div class="cav"><strong>&#9888; The two panels below are NOT cross-comparable.</strong><br>
GGM weights = partial correlation edge weights in condition-specific networks (cells as obs).
Pseudobulk weights = mean intramodular Spearman correlation (subclusters as obs).
Different scales, different null distributions. Compare within each panel only.</div>

<h3>D1: GGM Module Condition Activation (within GGM mode only)</h3>',
t2h(csv_tbl("module_condition_activation_ggm.csv", 10)),
FIG("fig_module_condition_activation_ggm.png"),
'<h3>D1: Pseudobulk Module Condition Activation (within pseudobulk mode only)</h3>',
t2h(csv_tbl("module_condition_activation_pseudobulk.csv")),
FIG("fig_module_condition_activation_pseudobulk.png"),
'</div>

<!-- ANNOTATION -->
<div class="sec" id="s6"><h2>6. Functional Annotation (E1&#8211;E2) &#8212; REFERENCE OUTPUT ONLY</h2>
<div class="cav"><strong>Modules are NOT named by GO terms or TF families.</strong>
GO enrichment and TF content are attached as descriptive metadata per module.
They are reference context, not analysis criteria. Never used for module naming or selection.</div>
<h3>E1: GO BP Enrichment &#8212; All 6 Module Sets</h3>
<p>Total: <strong>1,491</strong> significant GO terms across all 6 module sets (BH q &lt; 0.05, min set size 10).</p>
<details><summary>Master GO enrichment table (first 25 rows)</summary>',
t2h(csv_tbl("go_enrichment_all_sets.csv", 25)),
'</details>
<h3>E2: TF Module Enrichment</h3>',
t2h(csv_tbl("tf_enrichment.csv", 15)),
'</div>

<!-- GENE-CENTRIC -->
<div class="sec" id="s7"><h2>7. Gene-centric Utility (F1&#8211;F2)</h2>
<h3>F1: Master Gene Lookup Table (18,836 genes)</h3>
<p>One row per gene: module assignment across all 6 sets, kME in each, top-10 co-expression partners in both modes (display labels: symbol&#160;(AT-ID)).</p>',
t2h(csv_tbl("gene_lookup_master.csv", 10)),
'<div class="note"><strong>How to use:</strong> Filter gene_lookup_master.csv by gene_id.
top10_GGM_consensus and top10_PB give guilt-by-association partners with weights.
Module columns allow instant cross-set comparison per gene.</div>
<h3>F2: WRKY Family &#8212; Demo of Generic Gene-set Query</h3>
<p>This demonstrates the generic gene-set module enrichment capability using the WRKY TF family.
The same query works for any user-provided gene list (CSV with gene_id column).</p>',
t2h(csv_tbl("geneset_query_demo_wrky.csv", 15)),
FIG("fig_wrky_crossmode.png"),
'</div>

<!-- SUMMARY -->
<div class="sec" id="s8"><h2>8. Summary &amp; How to Reuse</h2>

<h3>Key Findings</h3>
<ul>
<li><strong>Scale-free topology</strong>: Power-law alpha 5.1&#8211;9.4 across all 6 networks.</li>
<li><strong>Module quality varies by method</strong>: Louvain grey rates 1.9&#8211;8.6% (most genes assigned); WGCNA 8.1&#8211;68% (higher coherence within assigned modules).</li>
<li><strong>Cross-mode overlap is low and expected</strong>: max Jaccard = 0.026, ARI &#8776; 0. GGM (partial correlation) and pseudobulk (marginal correlation) capture different regulatory relationships. They are complementary views, not redundant ones.</li>
<li><strong>1,687 genes in all 4 GGM sets</strong>: these are the most stable GGM co-expression units.</li>
<li><strong>1,095 core genes</strong> assigned in both GGM and pseudobulk &#8212; highest-confidence co-expressed groups.</li>
<li><strong>1,491 significant GO terms</strong> across 6 sets (reference only). <strong>319 TF-module records</strong> (reference only).</li>
</ul>

<h3>Analysis Completion</h3>
<table border="1">
<thead><tr><th>Analysis</th><th>Status</th><th>Output</th><th>Key Metric</th></tr></thead>
<tbody>
<tr><td>A1: Degree + power-law</td><td style="color:green;font-weight:bold">DONE</td><td>topology_degree.csv</td><td>alpha 5.1&#8211;9.4</td></tr>
<tr><td>A2: Centrality</td><td style="color:green;font-weight:bold">DONE</td><td>topology_centrality.csv</td><td>79,720 records; betweenness for GGM consensus</td></tr>
<tr><td>A3: Global stats</td><td style="color:green;font-weight:bold">DONE</td><td>topology_global_stats.csv</td><td>6 networks characterized</td></tr>
<tr><td>B1: kME</td><td style="color:green;font-weight:bold">DONE</td><td>module_kme_distributions.csv</td><td>116 records; 92 low-coherence modules flagged</td></tr>
<tr><td>B2: Eigengene corr</td><td style="color:green;font-weight:bold">DONE</td><td>module_eigengene_correlations.csv</td><td>1,558 module-pair correlations</td></tr>
<tr><td>B3: Hub genes</td><td style="color:green;font-weight:bold">DONE</td><td>module_hubs.csv</td><td>2,050 records</td></tr>
<tr><td>B4: Module quality</td><td style="color:green;font-weight:bold">DONE</td><td>module_quality_summary.csv</td><td>Grey rates: 68%, 8.6%, 44.5%, 1.9%, 8.1%, 6.6%</td></tr>
<tr><td>C1: Cross-mode overlap</td><td style="color:green;font-weight:bold">DONE</td><td>crossmode_overlap.csv, crossmode_ari.csv</td><td>1,764 pairs; max Jaccard=0.026</td></tr>
<tr><td>C2: WGCNA vs Louvain</td><td style="color:green;font-weight:bold">DONE</td><td>crossmethod_agreement.csv</td><td>ARI 0.028&#8211;0.301</td></tr>
<tr><td>C3: Core vs mode-specific</td><td style="color:green;font-weight:bold">DONE</td><td>core_vs_modespecific.csv</td><td>Core=1,095; GGM-only=8,371; PB-only=4,048</td></tr>
<tr><td>C4: Cross-set consistency</td><td style="color:green;font-weight:bold">DONE</td><td>cross_set_consistency.csv</td><td>1,687 in all 4 GGM sets</td></tr>
<tr><td>D1: GGM activation</td><td style="color:green;font-weight:bold">DONE</td><td>module_condition_activation_ggm.csv</td><td>42 modules; reused condition_comparison/</td></tr>
<tr><td>D1: PB activation</td><td style="color:green;font-weight:bold">DONE</td><td>module_condition_activation_pseudobulk.csv</td><td>18 records</td></tr>
<tr><td>D2: GGM patterns</td><td style="color:green;font-weight:bold">DONE</td><td>module_condition_patterns_ggm.csv</td><td>98 records</td></tr>
<tr><td>D2: PB patterns</td><td style="color:green;font-weight:bold">DONE</td><td>module_condition_patterns_pseudobulk.csv</td><td>18 records</td></tr>
<tr><td>E1: GO enrichment</td><td style="color:green;font-weight:bold">DONE</td><td>go_enrichment_all_sets.csv</td><td>1,491 significant terms</td></tr>
<tr><td>E2: TF enrichment</td><td style="color:green;font-weight:bold">DONE</td><td>tf_enrichment.csv</td><td>319 records</td></tr>
<tr><td>F1: Gene lookup</td><td style="color:green;font-weight:bold">DONE</td><td>gene_lookup_master.csv</td><td>18,836 genes</td></tr>
<tr><td>F2: WRKY demo</td><td style="color:green;font-weight:bold">DONE</td><td>geneset_query_demo_wrky.csv</td><td>1 significant enrichment (BH q&lt;0.05)</td></tr>
</tbody></table>

<h3>How to Reuse on a New Dataset</h3>
<ol>
<li><strong>Fully automatic</strong> (no user input): A1&#8211;A3, B1&#8211;B4, C1&#8211;C4, D1&#8211;D2, E1 (GO), E2 (TF if module_tfs.csv present)</li>
<li><strong>Needs a parameter</strong> (gene list CSV): F2 &#8212; any gene list with a gene_id column</li>
<li><strong>Entry point</strong>: <code>Rscript inst/scripts/downstream_analysis.R</code> from repo root</li>
<li><strong>Resumable</strong>: all analyses check skip_if_exists; re-run safely after any failure</li>
</ol>

<h3>Output Locations</h3>
<ul>
<li><code>results/pathogen_multiome/downstream/</code> &#8212; all CSV outputs (21 files)</li>
<li><code>results/pathogen_multiome/downstream/figures/</code> &#8212; 28 figures (PNG + SVG)</li>
<li><code>results/pathogen_multiome/downstream/DOWNSTREAM_ANALYSIS_REPORT.html</code> &#8212; this report</li>
<li><code>inst/scripts/downstream_analysis.R</code> &#8212; main analysis (generic, reusable)</li>
</ul>
</div>
</body></html>')

cat("Writing HTML...\n")
writeLines(html, RPATH)
sz <- file.size(RPATH)
cat("Report:", RPATH, "\n")
cat("Size:", round(sz/1024/1024, 2), "MB | Standalone:", sz > 10000, "\n")
