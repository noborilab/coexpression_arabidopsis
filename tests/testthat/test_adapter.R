library(testthat)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Create a minimal Seurat object in memory, save to tempfile, return the path.
# Caller is responsible for cleanup via on.exit(unlink(path)).
.make_seurat_rds <- function(counts_mat, meta_df) {
  skip_if_not_installed("SeuratObject")
  obj <- SeuratObject::CreateSeuratObject(counts = counts_mat, meta.data = meta_df)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(obj, tmp)
  tmp
}

# 4 AT-ID genes x n_cells cells; every cell expresses every gene (count = 2).
.at_id_counts <- function(n_cells = 20) {
  matrix(
    2L, nrow = 4, ncol = n_cells,
    dimnames = list(
      c("AT1G01010", "AT2G34567", "AT3G12345", "AT4G98765"),
      paste0("cell_", seq_len(n_cells))
    )
  )
}

# Standard 20-cell metadata: condition = Mock (10) / DC3000 (10).
.std_meta <- function(n_cells = 20) {
  data.frame(
    condition = rep(c("Mock", "DC3000"), each = n_cells / 2L),
    sample_id = paste0("S", rep(seq_len(5L), times = n_cells / 5L)),
    row.names = paste0("cell_", seq_len(n_cells)),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. AT-ID row names pass through without mapping
# ---------------------------------------------------------------------------
test_that("AT-ID row names pass through without mapping", {
  skip_if_not_installed("SeuratObject")

  tmp <- .make_seurat_rds(.at_id_counts(), .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  result <- load_seurat(tmp, "test", stratum_var = "condition", slot = "counts")

  expect_setequal(rownames(result$counts),
                  c("AT1G01010", "AT2G34567", "AT3G12345", "AT4G98765"))
  expect_true(all(is.na(result$gene_meta$gene_symbol)))
  expect_setequal(result$gene_meta$gene_id, rownames(result$counts))
})

# ---------------------------------------------------------------------------
# 2. Symbol row names + symbol_map → correct mapping, unmapped dropped + warned
# ---------------------------------------------------------------------------
test_that("symbol_map maps correctly and warns about unmapped genes", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    2L, nrow = 4, ncol = 20,
    dimnames = list(
      c("FLS2", "WRKY33", "NPR1", "UNMAPPED_GENE"),
      paste0("cell_", seq_len(20))
    )
  )
  sym_map <- data.frame(
    gene_symbol = c("FLS2", "WRKY33", "NPR1"),
    gene_id     = c("AT5G46330", "AT2G38470", "AT1G64280"),
    stringsAsFactors = FALSE
  )

  tmp <- .make_seurat_rds(counts, .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  expect_warning(
    result <- load_seurat(tmp, "test", stratum_var = "condition",
                          slot = "counts", symbol_map = sym_map),
    "1 genes dropped"
  )

  expect_equal(nrow(result$gene_meta), 3L)
  expect_setequal(result$gene_meta$gene_id,
                  c("AT5G46330", "AT2G38470", "AT1G64280"))
  expect_setequal(result$gene_meta$gene_symbol, c("FLS2", "WRKY33", "NPR1"))
  expect_setequal(rownames(result$counts),
                  c("AT5G46330", "AT2G38470", "AT1G64280"))
})

# ---------------------------------------------------------------------------
# 3. Symbol row names, no symbol_map → warning emitted, symbols used as gene_id
# ---------------------------------------------------------------------------
test_that("no symbol_map emits a warning and uses symbols as gene_id", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    2L, nrow = 3, ncol = 20,
    dimnames = list(c("FLS2", "WRKY33", "NPR1"),
                    paste0("cell_", seq_len(20)))
  )

  tmp <- .make_seurat_rds(counts, .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  expect_warning(
    result <- load_seurat(tmp, "test", stratum_var = "condition",
                          slot = "counts"),
    "no symbol_map provided"
  )

  expect_setequal(result$gene_meta$gene_id, c("FLS2", "WRKY33", "NPR1"))
  expect_setequal(rownames(result$counts), c("FLS2", "WRKY33", "NPR1"))
})

# ---------------------------------------------------------------------------
# 4. stratum_var not in metadata → informative stop()
# ---------------------------------------------------------------------------
test_that("missing stratum_var raises an informative error", {
  skip_if_not_installed("SeuratObject")

  tmp <- .make_seurat_rds(.at_id_counts(), .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    load_seurat(tmp, "test", stratum_var = "NONEXISTENT", slot = "counts"),
    "stratum_var 'NONEXISTENT' not found"
  )
})

# ---------------------------------------------------------------------------
# 5. min_cells filter retains the correct gene count
# ---------------------------------------------------------------------------
test_that("min_cells filter retains only sufficiently expressed genes", {
  skip_if_not_installed("SeuratObject")

  # 5 genes x 20 cells
  # Genes 1-3: expressed in all 20 cells
  # Genes 4-5: expressed in only 5 cells
  counts <- matrix(0L, nrow = 5L, ncol = 20L,
    dimnames = list(
      c("AT1G01010", "AT2G34567", "AT3G12345", "AT4G98765", "AT5G11111"),
      paste0("cell_", seq_len(20))
    ))
  counts[1:3, ]    <- 2L
  counts[4:5, 1:5] <- 2L

  tmp <- .make_seurat_rds(counts, .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  # min_cells = 10: only genes 1-3 pass (20 cells ≥ 10)
  res_strict <- load_seurat(tmp, "test", stratum_var = "condition",
                             slot = "counts", min_cells = 10L)
  expect_equal(nrow(res_strict$counts), 3L)

  # min_cells = 5: all 5 genes pass (genes 4-5 have exactly 5 expressing cells)
  res_loose <- load_seurat(tmp, "test", stratum_var = "condition",
                            slot = "counts", min_cells = 5L)
  expect_equal(nrow(res_loose$counts), 5L)
})

# ---------------------------------------------------------------------------
# 6. stratum_levels subsetting retains the correct cell count
# ---------------------------------------------------------------------------
test_that("stratum_levels subsets cells to the requested levels only", {
  skip_if_not_installed("SeuratObject")

  tmp <- .make_seurat_rds(.at_id_counts(), .std_meta())  # 10 Mock, 10 DC3000
  on.exit(unlink(tmp), add = TRUE)

  result <- load_seurat(tmp, "test", stratum_var = "condition",
                        stratum_levels = "Mock", slot = "counts")

  expect_equal(ncol(result$counts), 10L)
  expect_true(all(result$cell_meta$condition == "Mock"))
  expect_equal(result$stratum_spec$levels, "Mock")
})

# ---------------------------------------------------------------------------
# 7. Output is a named list with all five required InputBundle slots
# ---------------------------------------------------------------------------
test_that("output is a named list with the five required InputBundle slots", {
  skip_if_not_installed("SeuratObject")

  tmp <- .make_seurat_rds(.at_id_counts(), .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  result <- load_seurat(tmp, "pathogen_multiome", stratum_var = "condition",
                        slot = "counts")

  expect_type(result, "list")
  expect_named(result,
               c("counts", "cell_meta", "gene_meta", "stratum_spec", "dataset_id"),
               ignore.order = TRUE)
  expect_equal(result$dataset_id, "pathogen_multiome")
  expect_s3_class(result$cell_meta, "data.frame")
  expect_s3_class(result$gene_meta, "data.frame")
  expect_true("cell_id" %in% names(result$cell_meta))
  expect_named(result$stratum_spec, c("variable", "levels"))
  expect_equal(result$stratum_spec$variable, "condition")
  expect_true(all(c("gene_id", "gene_symbol") %in% names(result$gene_meta)))
})

# ---------------------------------------------------------------------------
# 8. counts rownames match gene_meta$gene_id in all cases
# ---------------------------------------------------------------------------
test_that("counts rownames are AT-IDs and match gene_meta$gene_id", {
  skip_if_not_installed("SeuratObject")

  tmp <- .make_seurat_rds(.at_id_counts(), .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  result <- load_seurat(tmp, "test", stratum_var = "condition", slot = "counts")

  expect_setequal(rownames(result$counts), result$gene_meta$gene_id)
  # All rownames should match the AT-ID pattern
  expect_true(all(grepl("^AT[1-5MC]G[0-9]{5}$", rownames(result$counts))))
})

test_that("counts rownames are symbols (not AT-IDs) when no mapping provided", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    2L, nrow = 2, ncol = 20,
    dimnames = list(c("FLS2", "BAK1"), paste0("cell_", seq_len(20)))
  )

  tmp <- .make_seurat_rds(counts, .std_meta())
  on.exit(unlink(tmp), add = TRUE)

  expect_warning(
    result <- load_seurat(tmp, "test", stratum_var = "condition",
                          slot = "counts"),
    "no symbol_map provided"
  )
  expect_setequal(rownames(result$counts), result$gene_meta$gene_id)
  expect_setequal(rownames(result$counts), c("FLS2", "BAK1"))
})
