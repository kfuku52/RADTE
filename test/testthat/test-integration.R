# Integration tests: run RADTE end-to-end and verify output

radte_script <- file.path(project_root, "radte.r")

test_that("RADTE NOTUNG mode produces expected output files", {
  sp_tree_file <- file.path(project_root, "data", "example_notung_01", "species_tree.nwk")
  gn_tree_file <- file.path(project_root, "data", "example_notung_01", "gene_tree.nwk.reconciled")
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(all(file.exists(c(sp_tree_file, gn_tree_file, parsable_file))),
              "NOTUNG example data not found")

  out_dir <- file.path(tempdir(), "radte_test_notung")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--gene_tree=", shQuote(gn_tree_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000",
    "--chronos_lambda=1",
    "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # Check expected output files exist
  expect_true(file.exists(file.path(out_dir, "radte_gene_tree_output.nwk")))
  expect_true(file.exists(file.path(out_dir, "radte_calibrated_nodes.txt")))
  expect_true(file.exists(file.path(out_dir, "radte_calibration_used.tsv")))
  expect_true(file.exists(file.path(out_dir, "radte_calibration_all.tsv")))
  expect_true(file.exists(file.path(out_dir, "radte_gene_tree.tsv")))
  expect_true(file.exists(file.path(out_dir, "radte_species_tree.tsv")))

  # Check output tree is a valid Newick tree
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(inherits(out_tree, "phylo"))
  expect_gt(length(out_tree$tip.label), 0)

  # Check calibrated nodes file
  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_true(cal_nodes[1] %in% c("RS", "S", "R", "allS"))
})

test_that("RADTE GeneRax mode produces expected output files", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  sp_tree_file <- file.path(project_root, "data", "example_generax_01", "species_tree.nwk")
  nhx_file <- file.path(project_root, "data", "example_generax_01", "gene_tree.nhx")
  skip_if_not(all(file.exists(c(sp_tree_file, nhx_file))),
              "GeneRax example data not found")

  out_dir <- file.path(tempdir(), "radte_test_generax")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--generax_nhx=", shQuote(nhx_file)),
    "--max_age=1000",
    "--chronos_lambda=1",
    "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # Check expected output files exist
  expect_true(file.exists(file.path(out_dir, "radte_gene_tree_output.nwk")))
  expect_true(file.exists(file.path(out_dir, "radte_calibrated_nodes.txt")))
  expect_true(file.exists(file.path(out_dir, "radte_calibration_used.tsv")))

  # Check output tree is a valid Newick tree
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(inherits(out_tree, "phylo"))
  expect_gt(length(out_tree$tip.label), 0)
})
