test_that("generated CLI exposes help and version without loading inputs", {
  help_output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(radte_path), "--help"),
    stdout = TRUE,
    stderr = TRUE
  )
  version_output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(radte_path), "--version"),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(help_output, "status"))
  expect_match(paste(help_output, collapse = "\n"), "Usage:")
  expect_null(attr(version_output, "status"))
  expect_identical(tail(version_output, 1L), radte_version)
})

test_that("atomic output keeps the previous file when a writer fails", {
  target <- tempfile()
  writeLines("before", target)
  on.exit(unlink(target))

  expect_error(
    atomic_write_file(target, function(tmp) stop("simulated writer failure")),
    "simulated writer failure"
  )
  expect_identical(readLines(target), "before")
})

test_that("SHA-256 and run manifest helpers produce auditable output", {
  input <- tempfile()
  manifest_file <- tempfile(fileext = ".tsv")
  writeLines("RADTE", input)
  on.exit(unlink(c(input, manifest_file)))

  digest <- compute_file_sha256(input)
  expect_match(digest, "^[0-9a-f]{64}$")

  write_run_manifest(c(version = radte_version, input_sha256 = digest), manifest_file)
  manifest <- read.delim(manifest_file, colClasses = "character")
  expect_equal(manifest$key, c("version", "input_sha256"))
  expect_equal(manifest$value[[2]], digest)
})

test_that("canonical clade matching is topology- and order-aware", {
  source_tree <- read.tree(text = "((A:1,B:1)ab:1,C:2)root;")
  target_tree <- read.tree(text = "(C:2,(B:1,A:1)x:1)y;")
  transferred <- transfer_node_labels(source_tree, target_tree)
  expect_equal(transferred$node.label, c("root", "ab"))

  different_tree <- read.tree(text = "((A:1,C:1)x:1,B:2)y;")
  different_tree$node.label[] <- "unchanged"
  transferred_different <- transfer_node_labels(source_tree, different_tree)
  expect_equal(transferred_different$node.label[[2]], "unchanged")
})
