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

  # Regression: output tree should be ultrametric
  expect_true(is.ultrametric(out_tree))

  # Regression: output tree should preserve tip labels from input
  in_tree <- read.tree(gn_tree_file)
  expect_setequal(out_tree$tip.label, in_tree$tip.label)

  # Regression: calibration_used.tsv should have valid structure
  cal_used <- read.delim(file.path(out_dir, "radte_calibration_used.tsv"))
  expect_true("node" %in% colnames(cal_used))
  expect_true("age.min" %in% colnames(cal_used))
  expect_true("age.max" %in% colnames(cal_used))
  expect_true(all(cal_used$age.min <= cal_used$age.max))

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

  # Regression: output tree should be ultrametric
  expect_true(is.ultrametric(out_tree))

  # Regression: calibration_used.tsv should have valid structure
  cal_used <- read.delim(file.path(out_dir, "radte_calibration_used.tsv"))
  expect_true("node" %in% colnames(cal_used))
  expect_true("age.min" %in% colnames(cal_used))
  expect_true("age.max" %in% colnames(cal_used))
  expect_true(all(cal_used$age.min <= cal_used$age.max))
})

# --- Regression: NOTUNG mode with no duplication events (allS mode) ---

test_that("RADTE NOTUNG mode handles gene tree without duplications (allS)", {
  # Create species tree
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((Homo_sapiens:10,Mus_musculus:10)n1:20,Danio_rerio:30)root;", sp_file)
  on.exit(unlink(sp_file))

  # Create gene tree with only speciation events (one gene per species)
  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((Homo_sapiens_g1:0.1,Mus_musculus_g1:0.1)n1:0.2,Danio_rerio_g1:0.3)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  # Parsable file with NO duplication events
  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t3\t5\t3\t0.1,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_test_allS")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
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

  # Should produce allS calibrated nodes
  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_equal(cal_nodes[1], "allS")

  # Output tree should exist and be valid
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(inherits(out_tree, "phylo"))
  expect_true(is.ultrametric(out_tree))
})
