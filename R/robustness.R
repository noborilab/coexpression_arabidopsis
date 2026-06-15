#' @title Cross-Context Robustness Layer (optional)
#'
#' @description
#' Computes cross-stratum and cross-dataset reproducibility scores for
#' co-expression edges. Enabled or disabled via `robustness.enabled` in config.
#'
#' **R_score method** (per-stratum fixed-evidence indicator aggregation):
#' - Per-stratum Fisher z-transform: `z_s = atanh(rho_s)`, `SE_s = 1/sqrt(n_s - 3)`
#' - Fixed-evidence indicator: `I_s = 1[z_s >= k * SE_s]` where `k ~ 1.64`
#'   (calibrated on positive controls; small strata require larger rho to qualify)
#' - Weighted aggregate: `R_score = sum(w_s * I_s) / sum(w_s)`,
#'   `w_s = sqrt(min(n_s, 30) - 3)`
#' - Null: analytic normal approximation to weighted Poisson-binomial
#'   (pi_s = P(Z >= k) under H0, same for all strata) → BH-FDR.
#'
#' **Benchmark note:** R_score did not outperform naive Spearman on GO co-functional
#' pair recovery (AUPRC ~0.20 for all methods). Value is interpretability and
#' context annotation, not global network improvement.
#'
#' @name robustness
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Canonical pair key: smaller AT-ID first, separated by TAB.
# Vectorised — never iterates row-by-row.
.canonical_key <- function(a, b) {
  swap  <- a > b
  key_a <- ifelse(swap, b, a)
  key_b <- ifelse(swap, a, b)
  paste(key_a, key_b, sep = "\t")
}

# Extract the sample size n_s from a NetworkResult$params list.
# GGM stores n_cells; pseudobulk will store n_pseudobulk.
.extract_n <- function(params) {
  n <- params[["n_cells"]] %||% params[["n_pseudobulk"]]
  if (is.null(n) || !is.numeric(n) || n < 4) {
    stop("NetworkResult$params must contain 'n_cells' or 'n_pseudobulk' >= 4.")
  }
  as.numeric(n)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------------------------------------------------------------------------

#' Compute cross-stratum robustness statistics
#'
#' @param network_list Named list of NetworkResult (one per stratum), from
#'   [estimate_pseudobulk()] or [estimate_singlecellggm()].
#' @param k Z-score multiplier for the fixed-evidence indicator I_s.
#'   Default 1.64 (one-sided 95th percentile). Larger k = more stringent
#'   per-stratum evidence required.
#' @param weight_cap Cap on the stratum weight: w_s = sqrt(min(n_s, weight_cap) - 3).
#'   Default 30. Prevents large-cell strata from dominating in GGM mode.
#' @param fdr_method FDR method passed to [p.adjust()]. Default "BH".
#' @return A RobustnessResult (named list; see `docs/OUTPUT_SCHEMA.md`).
#'   `$pair_scores` contains ALL tested pairs (FLAG-03: do not filter here).
#'   `$method_params` records all parameters used.
#' @export
compute_robustness <- function(network_list,
                               k          = 1.64,
                               weight_cap = 30,
                               fdr_method = "BH") {

  if (!is.list(network_list) || is.null(names(network_list)) ||
      any(nchar(names(network_list)) == 0)) {
    stop("network_list must be a named list of NetworkResult.")
  }

  strata <- names(network_list)
  S      <- length(strata)

  # ---- per-stratum metadata -----------------------------------------------

  n_s <- vapply(strata, function(s) .extract_n(network_list[[s]]$params),
                numeric(1))

  w_s     <- sqrt(pmax(0, pmin(n_s, weight_cap) - 3))
  sum_w   <- sum(w_s)

  if (sum_w == 0) stop("All stratum weights are zero (n_s <= 3 in all strata).")

  # Fixed null indicator probability: P(Z >= k) under H0 (analytic; same for
  # all strata because the SE_s cancels in the standardised threshold).
  # TODO: exact permutation-based pi_s per (n_s, bin) for better calibration.
  pi_s <- pnorm(k, lower.tail = FALSE)

  # ---- build union of all tested pairs ------------------------------------

  # One row per edge per stratum; immediately compute canonical key.
  pair_chunks <- lapply(strata, function(s) {
    et <- network_list[[s]]$edge_table
    if (is.null(et) || nrow(et) == 0L) return(NULL)
    data.frame(
      key    = .canonical_key(et$gene_id_A, et$gene_id_B),
      weight = et$weight,
      s      = s,
      stringsAsFactors = FALSE
    )
  })
  pair_all <- do.call(rbind, Filter(Negate(is.null), pair_chunks))

  if (is.null(pair_all) || nrow(pair_all) == 0L) {
    stop("No edges found across any stratum in network_list.")
  }

  unique_keys <- unique(pair_all$key)
  n_pairs     <- length(unique_keys)

  # Recover ordered (A, B) from key
  key_parts <- strsplit(unique_keys, "\t", fixed = TRUE)
  pair_A    <- vapply(key_parts, `[`, character(1), 1)
  pair_B    <- vapply(key_parts, `[`, character(1), 2)

  # ---- per-pair, per-stratum matrices -------------------------------------
  # rho_mat[i, s] = weight from edge_table (0 if pair absent in stratum s)
  # z_mat[i, s]   = atanh(rho_mat[i, s])
  # I_mat[i, s]   = 1 if z_mat[i,s] >= k * SE_s else 0

  rho_mat <- matrix(0.0, nrow = n_pairs, ncol = S)
  z_mat   <- matrix(0.0, nrow = n_pairs, ncol = S)
  I_mat   <- matrix(0L,  nrow = n_pairs, ncol = S)

  for (si in seq_along(strata)) {
    s         <- strata[si]
    sub_df    <- pair_chunks[[si]]
    if (is.null(sub_df)) next

    SE_s  <- 1.0 / sqrt(pmax(1, n_s[s] - 3))
    thresh <- k * SE_s

    lookup  <- setNames(sub_df$weight, sub_df$key)
    rho_s   <- lookup[unique_keys]   # NA for pairs absent this stratum
    rho_s[is.na(rho_s)] <- 0.0

    z_s <- atanh(rho_s)

    rho_mat[, si] <- rho_s
    z_mat[, si]   <- z_s
    I_mat[, si]   <- as.integer(z_s >= thresh)
  }

  # ---- aggregate statistics -----------------------------------------------

  R_score <- as.numeric(I_mat %*% w_s) / sum_w
  z_bar   <- as.numeric(z_mat %*% w_s) / sum_w

  # Between-stratum heterogeneity (DerSimonian-Laird estimator):
  #   Q     = sum_s w_s * (z_s - z_bar)^2
  #   tau2  = max(0, (Q - (S-1)) / (sum_w - sum(w_s^2) / sum_w))
  z_dev  <- z_mat - matrix(z_bar, nrow = n_pairs, ncol = S, byrow = FALSE)
  Q      <- as.numeric((z_dev^2) %*% w_s)
  denom  <- sum_w - sum(w_s^2) / sum_w
  tau2   <- pmax(0, (Q - (S - 1L)) / if (denom > 0) denom else 1)

  # p-value via normal approximation to weighted Poisson-binomial null:
  mu_null  <- sum(w_s * pi_s) / sum_w          # same for all pairs
  var_null <- sum(w_s^2 * pi_s * (1 - pi_s)) / sum_w^2

  z_null <- (R_score - mu_null) / sqrt(var_null)
  pval   <- pnorm(z_null, lower.tail = FALSE)
  qval   <- p.adjust(pval, method = fdr_method)

  # ---- assemble pair_scores -----------------------------------------------

  pair_scores <- data.frame(
    gene_id_A = pair_A,
    gene_id_B = pair_B,
    R_score   = R_score,
    z_bar     = z_bar,
    tau2      = tau2,
    pval      = pval,
    qval      = qval,
    stringsAsFactors = FALSE
  )

  for (si in seq_along(strata)) {
    pair_scores[[paste0("I_", strata[si])]] <- I_mat[, si]
  }

  pair_scores$star <- NA_real_   # filled later by annotate_star()

  list(
    pair_scores   = pair_scores,
    method_params = list(
      k             = k,
      weight_cap    = weight_cap,
      fdr_method    = fdr_method,
      n_strata      = S,
      stratum_names = strata
    )
  )
}

# ---------------------------------------------------------------------------

#' Cross-dataset replication annotation
#'
#' Adds a `star` column to `rob1$pair_scores`:
#' `TRUE` if the pair has `R_score >= threshold` in **both** datasets,
#' `FALSE` if the pair is in both datasets but below threshold in rob2,
#' `NA` if the pair is absent from rob2 (no replication data available).
#'
#' @param rob1 RobustnessResult from dataset 1 (will be annotated).
#' @param rob2 RobustnessResult from dataset 2 (replication dataset; not modified).
#' @param threshold R_score threshold for "robust in both". Default 0.5.
#' @return `rob1` with `$pair_scores$star` filled.
#' @export
annotate_star <- function(rob1, rob2, threshold = 0.5) {

  ps1 <- rob1$pair_scores
  ps2 <- rob2$pair_scores

  # Build rob2 lookup: canonical key → R_score (handle both orientations)
  key2_fwd <- .canonical_key(ps2$gene_id_A, ps2$gene_id_B)
  r2_lookup <- setNames(ps2$R_score, key2_fwd)

  # Look up canonical keys from rob1
  key1   <- .canonical_key(ps1$gene_id_A, ps1$gene_id_B)
  r2_val <- r2_lookup[key1]          # NA when pair absent from rob2

  present  <- !is.na(r2_val)
  star     <- rep(NA_real_, nrow(ps1))
  star[present] <- as.numeric(
    ps1$R_score[present] >= threshold & r2_val[present] >= threshold
  )

  rob1$pair_scores$star <- star
  rob1
}

# ---------------------------------------------------------------------------

#' Save RobustnessResult to disk
#'
#' Writes `pair_scores_full.csv` (ALL pairs, never filtered — FLAG-03) and
#' `robustness_result.rds`.
#'
#' @param rob RobustnessResult from [compute_robustness()].
#' @param outdir Output directory (created if it does not exist).
#' @return `outdir` invisibly.
#' @export
save_robustness <- function(rob, outdir) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  csv_path <- file.path(outdir, "pair_scores_full.csv")
  rds_path <- file.path(outdir, "robustness_result.rds")

  write.csv(rob$pair_scores, csv_path, row.names = FALSE)
  saveRDS(rob, rds_path)

  message("Saved: ", nrow(rob$pair_scores), " pairs -> ", csv_path)
  message("Saved: RDS -> ", rds_path)

  invisible(outdir)
}

# ---------------------------------------------------------------------------

#' Characterize each gene pair by its cross-condition activity pattern
#'
#' Produces BOTH a discrete binary pattern label and the continuous
#' per-condition effect sizes, so pairs can be grouped by pattern OR clustered
#' on the continuous matrix. The discrete pattern is derived mechanically from
#' the I_s indicators already present in `rob$pair_scores`.
#'
#' **Pattern label behaviour:**
#' - When `pattern_labels = NULL` (default): all patterns receive a generic
#'   label of the form `"pattern_<bits>"` (e.g. `"pattern_1111"`,
#'   `"pattern_0001"`). No biology-specific strings are ever embedded in the
#'   library function.
#' - When `pattern_labels` is supplied: each pattern code is looked up in the
#'   named vector; any pattern not present in the lookup gets the generic
#'   `"pattern_<bits>"` label and a warning is emitted. Caller is responsible
#'   for supplying biologically meaningful labels appropriate to their dataset.
#'
#' **Pathogen multiome usage** (in the runner scripts):
#' ```r
#' PATTERN_LABELS <- c(
#'   "0000" = "none",          "1111" = "constitutive_all",
#'   "1000" = "single_Mock",   "0100" = "single_DC3000",
#'   "0010" = "single_AvrRpt2","0001" = "single_AvrRpm1",
#'   "0111" = "pan_pathogen",  "0011" = "ETI_shared"
#' )
#' cp <- characterize_condition_pattern(rob, network_list,
#'         condition_order = CONDITIONS, pattern_labels = PATTERN_LABELS)
#' ```
#'
#' **Specificity index:** `(w_max - w_mean_of_others) / (w_max + epsilon)`,
#' where `w_mean_of_others = (sum(w) - w_max) / (S - 1)` and `epsilon = 1e-6`.
#' Approaches 1 when signal is concentrated in one condition; approaches 0
#' when uniform across all conditions.
#'
#' @param rob RobustnessResult from [compute_robustness()]. `$pair_scores` must
#'   contain `I_<condition>` columns for each element of `condition_order`.
#' @param network_list Named list of NetworkResult (one per condition). Names
#'   must include each element of `condition_order`. Conditions absent from the
#'   list are treated as weight = 0 for all pairs.
#' @param condition_order Character vector fixing the bit order for the pattern
#'   string. Default `c("Mock","DC3000","AvrRpt2","AvrRpm1")`.
#' @param pattern_labels Optional named character vector mapping bit-pattern
#'   codes (e.g. `"0111"`) to human-readable labels. `NULL` (default) produces
#'   generic `"pattern_<bits>"` labels for all patterns. Any pattern not found
#'   in the lookup gets the generic label and triggers a warning.
#' @return `data.frame`, one row per pair in `rob$pair_scores`, with columns:
#'   `gene_id_A`, `gene_id_B`;
#'   `I_<condition>` (integer 0/1, from rob);
#'   `pattern` (S-character bit string in condition_order);
#'   `pattern_label` (category string; see `pattern_labels` parameter);
#'   `n_conditions_active` (integer 0–S);
#'   `w_<condition>` (pcor weight from per-condition network, or 0 if absent);
#'   `w_max`, `w_min`, `w_range`, `w_mean`, `specificity_index`.
#' @export
characterize_condition_pattern <- function(rob, network_list,
  condition_order = c("Mock", "DC3000", "AvrRpt2", "AvrRpm1"),
  pattern_labels  = NULL) {

  ps <- rob$pair_scores
  S  <- length(condition_order)

  i_cols <- paste0("I_", condition_order)
  w_cols <- paste0("w_", condition_order)

  # ---- continuous weights: look up per-condition edge weight ----------------

  pair_key <- .canonical_key(ps$gene_id_A, ps$gene_id_B)
  w_mat    <- matrix(0.0, nrow = nrow(ps), ncol = S,
                     dimnames = list(NULL, w_cols))

  for (ci in seq_len(S)) {
    cond <- condition_order[ci]
    nr   <- network_list[[cond]]
    if (is.null(nr)) next
    et <- nr$edge_table
    if (is.null(et) || nrow(et) == 0L) next
    ck       <- .canonical_key(et$gene_id_A, et$gene_id_B)
    w_lookup <- setNames(et$weight, ck)
    matched  <- w_lookup[pair_key]
    matched[is.na(matched)] <- 0.0
    w_mat[, ci] <- as.numeric(matched)
  }

  # ---- discrete pattern from I_s indicators in rob$pair_scores --------------

  I_mat <- matrix(NA_integer_, nrow = nrow(ps), ncol = S)
  for (ci in seq_len(S)) {
    col <- i_cols[ci]
    if (col %in% names(ps)) I_mat[, ci] <- as.integer(ps[[col]])
  }

  pattern <- apply(I_mat, 1L, paste, collapse = "")

  # Build label map from unique patterns.
  # When pattern_labels is supplied: use lookup; generic "pattern_<bits>" for unmapped.
  # When pattern_labels is NULL:     use generic "pattern_<bits>" for all patterns.
  unique_pats <- unique(pattern)
  if (!is.null(pattern_labels)) {
    unmapped <- setdiff(unique_pats, names(pattern_labels))
    if (length(unmapped) > 0L)
      warning("characterize_condition_pattern: ", length(unmapped),
              " pattern(s) not in pattern_labels; using generic 'pattern_<bits>': ",
              paste(unmapped, collapse = ", "))
    label_map <- setNames(
      vapply(unique_pats, function(p) {
        if (p %in% names(pattern_labels)) pattern_labels[[p]]
        else paste0("pattern_", p)
      }, character(1L)),
      unique_pats
    )
  } else {
    label_map <- setNames(paste0("pattern_", unique_pats), unique_pats)
  }
  pattern_label <- unname(label_map[pattern])

  n_conditions_active <- rowSums(I_mat, na.rm = TRUE)

  # ---- continuous statistics -----------------------------------------------

  w_max   <- apply(w_mat, 1L, max)
  w_min   <- apply(w_mat, 1L, min)
  w_range <- w_max - w_min
  w_mean  <- rowMeans(w_mat)

  epsilon       <- 1e-6
  w_mean_others <- (rowSums(w_mat) - w_max) / max(S - 1L, 1L)
  specificity_index <- (w_max - w_mean_others) / (w_max + epsilon)

  # ---- assemble result ------------------------------------------------------

  out <- data.frame(gene_id_A = ps$gene_id_A,
                    gene_id_B = ps$gene_id_B,
                    stringsAsFactors = FALSE)
  for (ci in seq_len(S)) out[[i_cols[ci]]] <- I_mat[, ci]
  out$pattern             <- pattern
  out$pattern_label       <- pattern_label
  out$n_conditions_active <- as.integer(n_conditions_active)
  for (ci in seq_len(S)) out[[w_cols[ci]]] <- w_mat[, ci]
  out$w_max             <- w_max
  out$w_min             <- w_min
  out$w_range           <- w_range
  out$w_mean            <- w_mean
  out$specificity_index <- specificity_index

  out
}
