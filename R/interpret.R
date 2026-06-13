#' @title Module Construction and Biological Interpretation
#'
#' @description
#' Shared output/interpretation layer used by **both** estimation modes.
#'
#' Covers:
#' - WGCNA-style signed-network module construction (soft power, merge threshold)
#' - Granularity sweep / hierarchy: coarse top-level modules with nested sub-modules
#' - Cross-context module preservation (WGCNA `modulePreservation` Zsummary,
#'   or lightweight mean intramodular |cor| z-score fallback when full
#'   `modulePreservation` times out on large matrices)
#' - Hub genes (kME)
#' - TF intersection (regulator hints per module)
#' - GO BP enrichment
#' - Curated-set anchor enrichment (fold-enrichment; WGCNA beat hclust-on-organ-means
#'   50× vs 27× on curated anchors in benchmarks)
#' - Plain-language per-module labels
#'
#' @name interpret
NULL

#' Build co-expression modules from a network
#'
#' @param network Gene × gene correlation or partial-correlation matrix.
#' @param soft_power Soft-thresholding power for signed network adjacency.
#'   If `NULL`, estimated automatically by scale-free topology fit.
#' @param merge_threshold Module merge height threshold (0–1). Controls
#'   granularity; higher values produce fewer, coarser modules.
#'
#' @return Named list:
#'   - `modules`: named integer vector of gene → module assignments
#'   - `kME`: data.frame of module eigengene correlations (genes × modules)
#'   - `dendro`: hierarchical clustering dendrogram object
#' @export
build_modules <- function(network, soft_power = NULL, merge_threshold = 0.25) {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}

#' Compute module preservation across contexts
#'
#' @param modules Output of [build_modules()].
#' @param reference_network The network on which modules were originally built.
#' @param test_networks Named list of networks to test preservation against.
#' @param method One of `"modulePreservation"` (WGCNA; accurate but slow —
#'   cap or subsample for matrices > ~1000 samples × 10k genes) or
#'   `"intramodular_cor"` (lightweight mean-intramodular-|cor| z-score fallback).
#'   Default `"intramodular_cor"`.
#'
#' @return `data.frame` of per-module preservation statistics (Zsummary or
#'   z-score proxy, p-value, effect size).
#' @export
compute_preservation <- function(modules,
                                 reference_network,
                                 test_networks,
                                 method = "intramodular_cor") {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}

#' Annotate modules with biological interpretation
#'
#' @param modules Output of [build_modules()].
#' @param go_annotations `data.frame` of GO BP annotations with columns
#'   `gene_id` and `go_term`.
#' @param tf_list Optional character vector of transcription factor gene IDs.
#'   (Lab TF list: `Athaliana_motifs_metadata.tsv`, 673 TFs, column `motif_id`.)
#' @param curated_sets Optional named list of curated gene sets for
#'   anchor fold-enrichment scoring.
#'
#' @return `modules` object augmented with GO enrichment results, TF
#'   intersection counts, and curated-set fold-enrichment per module.
#' @export
annotate_modules <- function(modules,
                              go_annotations,
                              tf_list     = NULL,
                              curated_sets = NULL) {
  # TODO (Phase 2): implement
  stop("Not implemented yet — scaffold stub only")
}
