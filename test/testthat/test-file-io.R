# --- read_notung_parsable ---

test_that("read_notung_parsable reads NOTUNG parsable file correctly", {
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(file.exists(parsable_file), "NOTUNG example data not found")

  result <- read_notung_parsable(parsable_file, mode = "D")
  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result), c("event", "gn_node", "lower_sp_node", "upper_sp_node"))
  # The example has 10 duplication events (lines starting with #D, excluding header)
  expect_equal(nrow(result), 10)
  expect_true(all(result$event == "D"))
})

test_that("read_notung_parsable returns empty data frame for unsupported mode", {
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(file.exists(parsable_file), "NOTUNG example data not found")

  result <- read_notung_parsable(parsable_file, mode = "X")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(colnames(result), c("event", "gn_node", "lower_sp_node", "upper_sp_node"))
})

# --- read_generax_nhx ---

test_that("read_generax_nhx reads GeneRax NHX file correctly", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  nhx_file <- file.path(project_root, "data", "example_generax_01", "gene_tree.nhx")
  skip_if_not(file.exists(nhx_file), "GeneRax example data not found")

  # read_generax_nhx writes a temp file in the working directory, use tempdir
  old_wd <- getwd()
  setwd(tempdir())
  on.exit(setwd(old_wd))

  result <- read_generax_nhx(nhx_file)
  expect_true(inherits(result, "treedata"))
  phylo <- result@phylo
  expect_true(inherits(phylo, "phylo"))
  expect_gt(length(phylo$tip.label), 0)
})
