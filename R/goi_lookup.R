#' @title Gene-of-Interest (GOI) Lookup Resource
#'
#' @description
#' Builds the primary collaborator-facing deliverable: a per-gene lookup table
#' for a user-supplied GOI list. For each GOI the table reports module
#' assignment, module membership (kME), hub flag, cross-context preservation,
#' and top co-expressed partners (FLAG-07: both wide and long formats).
#'
#' @name goi_lookup
NULL

#' Build the GOI lookup table
#'
#' For each gene in `goi_ids`, looks up its module, kME, hub status,
#' preservation score, and top co-expressed partners. Genes not in the network
#' receive a row with `NA` module and `notes = "not in network"`.
#'
#' @param mod_input ModuleInput from interpret.R functions.
#' @param rob RobustnessResult from [compute_robustness()] (for partner
#'   R_scores and z_bar weights).
#' @param goi_ids Character vector of AT-IDs to look up.
#' @param top_n Number of top co-expressed partners to report. Default `20`.
#' @return Named list with two elements:
#'   - `$wide`: one row per GOI.
#'   - `$long`: one row per GOI × partner pair (FLAG-07).
#' @export
build_goi_table <- function(mod_input, rob, goi_ids, top_n = 20L) {
  gene_mod  <- mod_input$gene_module
  mod_meta  <- mod_input$module_meta
  hub_genes <- mod_input$hub_genes
  ps        <- rob$pair_scores

  # Hub lookup: set of (gene_id, module_id) strings for O(1) lookup
  hub_key <- paste(hub_genes$gene_id, hub_genes$module_id)

  wide_list <- vector("list", length(goi_ids))
  long_list <- vector("list", length(goi_ids))

  for (i in seq_along(goi_ids)) {
    goi    <- goi_ids[i]
    gm_row <- gene_mod[gene_mod$gene_id == goi, , drop = FALSE]

    if (nrow(gm_row) == 0L) {
      # Gene not in the network
      wide_list[[i]] <- data.frame(
        gene_id                    = goi,
        gene_symbol                = NA_character_,
        module                     = NA_integer_,
        kME                        = NA_real_,
        hub_flag                   = FALSE,
        zsummary                   = NA_real_,
        preservation_method        = NA_character_,
        top_N_coexpressed_partners = NA_character_,
        notes                      = "not in network",
        stringsAsFactors = FALSE
      )
      long_list[[i]] <- NULL
      next
    }

    top_mod  <- gm_row$top_module[1L]
    kme_val  <- gm_row$kME[1L]
    hub_flag <- paste(goi, top_mod) %in% hub_key
    note     <- if (!is.na(top_mod) && top_mod == 0L) "unassigned" else ""

    # Module-level preservation
    if (!is.na(top_mod) && top_mod > 0L) {
      meta_row       <- mod_meta[mod_meta$module_id == top_mod, , drop = FALSE]
      zsummary_val   <- if (nrow(meta_row) > 0L) meta_row$zsummary[1L] else NA_real_
      pres_meth_val  <- if (nrow(meta_row) > 0L) meta_row$preservation_method[1L] else NA_character_
    } else {
      zsummary_val  <- NA_real_
      pres_meth_val <- NA_character_
    }

    # Top N co-expressed partners from RobustnessResult
    is_partner   <- (ps$gene_id_A == goi | ps$gene_id_B == goi)
    partner_rows <- ps[is_partner, , drop = FALSE]

    if (nrow(partner_rows) > 0L) {
      partner_rows$partner_id <- ifelse(
        partner_rows$gene_id_A == goi,
        partner_rows$gene_id_B,
        partner_rows$gene_id_A
      )
      partner_rows <- partner_rows[order(partner_rows$R_score, decreasing = TRUE), ]
      partner_rows <- head(partner_rows, top_n)
      top_n_str    <- paste(
        sprintf("%s (R=%.2f)", partner_rows$partner_id, partner_rows$R_score),
        collapse = "; "
      )
    } else {
      partner_rows <- partner_rows[0L, , drop = FALSE]  # empty, same structure
      top_n_str    <- ""
    }

    wide_list[[i]] <- data.frame(
      gene_id                    = goi,
      gene_symbol                = NA_character_,
      module                     = as.integer(top_mod),
      kME                        = kme_val,
      hub_flag                   = hub_flag,
      zsummary                   = zsummary_val,
      preservation_method        = pres_meth_val,
      top_N_coexpressed_partners = top_n_str,
      notes                      = note,
      stringsAsFactors = FALSE
    )

    if (nrow(partner_rows) > 0L) {
      long_list[[i]] <- data.frame(
        gene_id        = goi,
        gene_symbol    = NA_character_,
        partner_id     = partner_rows$partner_id,
        partner_symbol = NA_character_,
        R_score        = partner_rows$R_score,
        weight         = tanh(partner_rows$z_bar),
        rank           = seq_len(nrow(partner_rows)),
        notes          = note,
        stringsAsFactors = FALSE
      )
    } else {
      long_list[[i]] <- NULL
    }
  }

  wide_df <- do.call(rbind, wide_list)
  long_df <- do.call(rbind, Filter(Negate(is.null), long_list))

  if (is.null(long_df)) {
    long_df <- data.frame(
      gene_id = character(), gene_symbol = character(),
      partner_id = character(), partner_symbol = character(),
      R_score = numeric(), weight = numeric(), rank = integer(),
      notes = character(), stringsAsFactors = FALSE
    )
  }

  list(wide = wide_df, long = long_df)
}

#' Save GOI lookup table to disk
#'
#' Saves both the wide and long lookup tables as CSV files.
#'
#' @param goi_result Output of [build_goi_table()].
#' @param outdir Output directory (created if it does not exist).
#' @param list_name File name prefix, e.g. `"CZL_receptors"`.
#' @return `outdir` invisibly.
#' @export
save_goi_table <- function(goi_result, outdir, list_name) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  wide_path <- file.path(outdir, paste0(list_name, "_lookup_wide.csv"))
  long_path <- file.path(outdir, paste0(list_name, "_lookup_long.csv"))

  write.csv(goi_result$wide, wide_path, row.names = FALSE)
  write.csv(goi_result$long, long_path, row.names = FALSE)

  message("Saved: ", wide_path)
  message("Saved: ", long_path)

  invisible(outdir)
}
