library(testthat)

# ---------------------------------------------------------------------------
# Helper: write a minimal GFF3 to a temp file
# ---------------------------------------------------------------------------

.write_gff3 <- function(lines, gzip = FALSE) {
  path <- tempfile(fileext = if (gzip) ".gff3.gz" else ".gff3")
  if (gzip) {
    con <- gzfile(path, "w")
    writeLines(lines, con)
    close(con)
  } else {
    writeLines(lines, path)
  }
  path
}

GFF3_LINES <- c(
  "##gff-version 3",
  "# comment line",
  "1\taraport11\tgene\t3631\t5899\t.\t+\t.\tID=gene:AT1G01010;Name=NAC001;gene_id=AT1G01010",
  "1\taraport11\tmRNA\t3631\t5899\t.\t+\t.\tID=transcript:AT1G01010.1;Parent=gene:AT1G01010",
  "1\taraport11\tgene\t6788\t9130\t.\t-\t.\tID=gene:AT1G01020;Name=ARV1;gene_id=AT1G01020",
  "1\taraport11\tncRNA_gene\t11101\t11372\t.\t+\t.\tID=gene:AT1G03987;gene_id=AT1G03987",
  "1\taraport11\tgene\t11649\t13714\t.\t-\t.\tID=gene:AT1G01030;Name=NGA3;gene_id=AT1G01030",
  "1\taraport11\tgene\t23121\t31227\t.\t+\t.\tID=gene:AT1G01040;gene_id=AT1G01040"
)

# ---------------------------------------------------------------------------
# build_symbol_map tests
# ---------------------------------------------------------------------------

test_that("build_symbol_map returns correct columns", {
  path <- .write_gff3(GFF3_LINES)
  on.exit(unlink(path))
  sm <- build_symbol_map(path)
  expect_s3_class(sm, "data.frame")
  expect_named(sm, c("gene_id", "gene_symbol"))
})

test_that("build_symbol_map extracts gene_id and symbol correctly", {
  path <- .write_gff3(GFF3_LINES)
  on.exit(unlink(path))
  sm <- build_symbol_map(path)

  # Should have 4 gene features (not mRNA, not ncRNA_gene)
  expect_equal(nrow(sm), 4L)

  # Check gene IDs present
  expect_true("AT1G01010" %in% sm$gene_id)
  expect_true("AT1G01020" %in% sm$gene_id)
  expect_true("AT1G01030" %in% sm$gene_id)
  expect_true("AT1G01040" %in% sm$gene_id)

  # Check symbols
  expect_equal(sm$gene_symbol[sm$gene_id == "AT1G01010"], "NAC001")
  expect_equal(sm$gene_symbol[sm$gene_id == "AT1G01020"], "ARV1")
  expect_equal(sm$gene_symbol[sm$gene_id == "AT1G01030"], "NGA3")
})

test_that("build_symbol_map returns NA for genes without Name attribute", {
  path <- .write_gff3(GFF3_LINES)
  on.exit(unlink(path))
  sm <- build_symbol_map(path)
  expect_true(is.na(sm$gene_symbol[sm$gene_id == "AT1G01040"]))
})

test_that("build_symbol_map excludes non-gene features", {
  path <- .write_gff3(GFF3_LINES)
  on.exit(unlink(path))
  sm <- build_symbol_map(path)
  # mRNA and ncRNA_gene should not appear
  expect_false("AT1G01010.1" %in% sm$gene_id)
  expect_false("AT1G03987"   %in% sm$gene_id)
})

test_that("build_symbol_map errors on missing file", {
  expect_error(build_symbol_map("/nonexistent/path.gff3"), "not found")
})

test_that("build_symbol_map errors when no gene features found", {
  path <- .write_gff3(c("##gff-version 3",
                         "1\taraport11\tmRNA\t1\t100\t.\t+\t.\tID=t1"))
  on.exit(unlink(path))
  expect_error(build_symbol_map(path), "No 'gene' features")
})

test_that("build_symbol_map handles gzipped input", {
  path <- .write_gff3(GFF3_LINES, gzip = TRUE)
  on.exit(unlink(path))
  sm <- build_symbol_map(path)
  expect_equal(nrow(sm), 4L)
  expect_equal(sm$gene_symbol[sm$gene_id == "AT1G01010"], "NAC001")
})

# ---------------------------------------------------------------------------
# load_network_results / save_network_results round-trip
# ---------------------------------------------------------------------------

.make_network_list <- function(n = 3L) {
  ids <- paste0("AT1G", formatC(seq_len(n) * 10L, width = 5L, flag = "0"))
  pairs <- combn(n, 2)
  et <- data.frame(
    gene_id_A = ids[pairs[1, ]],
    gene_id_B = ids[pairs[2, ]],
    weight    = runif(ncol(pairs)),
    stringsAsFactors = FALSE
  )
  list(
    MockA = list(edge_table = et, gene_ids = ids, stratum_id = "MockA",
                 mode = "singlecellggm", params = list(n_cells = 100L),
                 timestamp = Sys.time())
  )
}

test_that("save and load network results round-trips correctly", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  nl <- .make_network_list()
  save_network_results(nl, tmpdir)
  nl2 <- load_network_results(tmpdir, strata = "MockA")

  expect_named(nl2, "MockA")
  expect_equal(nrow(nl2$MockA$edge_table), nrow(nl$MockA$edge_table))
  expect_equal(sort(colnames(nl2$MockA$edge_table)),
               sort(colnames(nl$MockA$edge_table)))
})

test_that("load_network_results errors on missing output_dir", {
  expect_error(load_network_results("/nonexistent/dir"), "does not exist")
})

test_that("save_network_results errors on unnamed list", {
  expect_error(save_network_results(list(1, 2), tempdir()),
               "named list")
})
