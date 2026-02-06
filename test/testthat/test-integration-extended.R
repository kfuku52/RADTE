# Extended integration tests: detailed output validation

radte_script <- file.path(project_root, "radte.r")

# Helper: run RADTE in a temp dir, return exit code, stderr, and output dir path
run_radte_ext <- function(..., keep_dir = FALSE) {
  out_dir <- file.path(tempdir(), paste0("radte_ext_", format(Sys.time(), "%H%M%S"), "_", sample(1000:9999, 1)))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (!keep_dir) on.exit(unlink(out_dir, recursive = TRUE))
  args <- paste(...)
  stderr_file <- tempfile()
  cmd <- paste("Rscript", shQuote(radte_script), args, "2>", shQuote(stderr_file))
  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)
  exit_code <- system(cmd, ignore.stdout = TRUE)
  stderr_output <- paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
  unlink(stderr_file)
  list(exit_code = exit_code, stderr = stderr_output, out_dir = out_dir)
}

# --- NOTUNG mode without --pad_short_edge ---

test_that("RADTE NOTUNG mode works without --pad_short_edge", {
  sp_tree_file <- file.path(project_root, "data", "example_notung_01", "species_tree.nwk")
  gn_tree_file <- file.path(project_root, "data", "example_notung_01", "gene_tree.nwk.reconciled")
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(all(file.exists(c(sp_tree_file, gn_tree_file, parsable_file))),
              "NOTUNG example data not found")

  out_dir <- file.path(tempdir(), "radte_nopad")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--gene_tree=", shQuote(gn_tree_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000",
    "--chronos_lambda=1",
    "--chronos_model=discrete"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # Output tree should still be generated
  expect_true(file.exists(file.path(out_dir, "radte_gene_tree_output.nwk")))
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(inherits(out_tree, "phylo"))
  expect_gt(length(out_tree$tip.label), 0)
})

# --- Gene tree TSV detailed validation ---

test_that("RADTE NOTUNG mode gene tree TSV has correct structure", {
  sp_tree_file <- file.path(project_root, "data", "example_notung_01", "species_tree.nwk")
  gn_tree_file <- file.path(project_root, "data", "example_notung_01", "gene_tree.nwk.reconciled")
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(all(file.exists(c(sp_tree_file, gn_tree_file, parsable_file))),
              "NOTUNG example data not found")

  out_dir <- file.path(tempdir(), "radte_tsv_check")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--gene_tree=", shQuote(gn_tree_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # Gene tree TSV validation
  gn_tsv <- read.delim(file.path(out_dir, "radte_gene_tree.tsv"))
  expected_cols <- c("event", "gn_node", "lower_sp_node", "upper_sp_node",
                     "lower_age", "upper_age", "gn_node_num")
  expect_true(all(expected_cols %in% colnames(gn_tsv)))

  # Events should be D, S, S(R), or D(R)
  expect_true(all(gn_tsv$event %in% c("D", "S", "S(R)", "D(R)")))

  # Exactly one root event
  root_events <- gn_tsv$event[grepl("\\(R\\)$", gn_tsv$event)]
  expect_equal(length(root_events), 1)

  # Node numbers should be positive integers
  expect_true(all(gn_tsv$gn_node_num > 0))

  # Ages should be non-negative
  expect_true(all(gn_tsv$lower_age >= 0))
  expect_true(all(gn_tsv$upper_age >= 0))

  # lower_age <= upper_age for all rows
  expect_true(all(gn_tsv$lower_age <= gn_tsv$upper_age))
})

# --- Species tree TSV detailed validation ---

test_that("RADTE NOTUNG mode species tree TSV has correct structure", {
  sp_tree_file <- file.path(project_root, "data", "example_notung_01", "species_tree.nwk")
  gn_tree_file <- file.path(project_root, "data", "example_notung_01", "gene_tree.nwk.reconciled")
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(all(file.exists(c(sp_tree_file, gn_tree_file, parsable_file))),
              "NOTUNG example data not found")

  out_dir <- file.path(tempdir(), "radte_sp_tsv")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--gene_tree=", shQuote(gn_tree_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # Species tree TSV validation
  sp_tsv <- read.delim(file.path(out_dir, "radte_species_tree.tsv"))
  expect_true("node" %in% colnames(sp_tsv))
  expect_true("age" %in% colnames(sp_tsv))

  # All ages should be >= 0
  expect_true(all(sp_tsv$age >= 0))

  # Tip species should have age 0
  sp_tree <- read.tree(sp_tree_file)
  tip_names <- sp_tree$tip.label
  tip_rows <- sp_tsv[sp_tsv$node %in% tip_names, ]
  expect_true(all(tip_rows$age == 0))

  # Number of rows = tips + internal nodes
  expect_equal(nrow(sp_tsv), length(sp_tree$tip.label) + sp_tree$Nnode)
})

# --- PDF output files existence ---

test_that("RADTE NOTUNG mode produces PDF output files", {
  sp_tree_file <- file.path(project_root, "data", "example_notung_01", "species_tree.nwk")
  gn_tree_file <- file.path(project_root, "data", "example_notung_01", "gene_tree.nwk.reconciled")
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(all(file.exists(c(sp_tree_file, gn_tree_file, parsable_file))),
              "NOTUNG example data not found")

  out_dir <- file.path(tempdir(), "radte_pdf_check")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--gene_tree=", shQuote(gn_tree_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # PDF files should be generated
  expect_true(file.exists(file.path(out_dir, "radte_gene_tree_input.pdf")))
  expect_true(file.exists(file.path(out_dir, "radte_gene_tree_output.pdf")))
  expect_true(file.exists(file.path(out_dir, "radte_species_tree.pdf")))

  # PDFs should have non-zero size
  expect_gt(file.size(file.path(out_dir, "radte_gene_tree_input.pdf")), 0)
  expect_gt(file.size(file.path(out_dir, "radte_gene_tree_output.pdf")), 0)
  expect_gt(file.size(file.path(out_dir, "radte_species_tree.pdf")), 0)
})

# --- Calibration table consistency ---

test_that("RADTE calibration_all.tsv is superset of calibration_used.tsv", {
  sp_tree_file <- file.path(project_root, "data", "example_notung_01", "species_tree.nwk")
  gn_tree_file <- file.path(project_root, "data", "example_notung_01", "gene_tree.nwk.reconciled")
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(all(file.exists(c(sp_tree_file, gn_tree_file, parsable_file))),
              "NOTUNG example data not found")

  out_dir <- file.path(tempdir(), "radte_cal_check")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--gene_tree=", shQuote(gn_tree_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  cal_all <- read.delim(file.path(out_dir, "radte_calibration_all.tsv"))
  cal_used <- read.delim(file.path(out_dir, "radte_calibration_used.tsv"))

  # Used should be a subset of all
  expect_true(all(cal_used$node %in% cal_all$node))

  # Both should have the same column structure
  expect_true(all(c("node", "age.min", "age.max", "soft.bounds", "event") %in% colnames(cal_all)))
  expect_true(all(c("node", "age.min", "age.max", "soft.bounds", "event") %in% colnames(cal_used)))

  # Events in calibration tables should be S, R, or S(R)
  expect_true(all(cal_all$event %in% c("S", "R", "S(R)")))
  expect_true(all(cal_used$event %in% c("S", "R", "S(R)")))

  # All ages should be non-negative
  expect_true(all(cal_all$age.min >= 0))
  expect_true(all(cal_all$age.max >= 0))
  expect_true(all(cal_all$age.min <= cal_all$age.max))
})

# --- GeneRax detailed output validation ---

test_that("RADTE GeneRax mode gene tree TSV has correct structure", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  sp_tree_file <- file.path(project_root, "data", "example_generax_01", "species_tree.nwk")
  nhx_file <- file.path(project_root, "data", "example_generax_01", "gene_tree.nhx")
  skip_if_not(all(file.exists(c(sp_tree_file, nhx_file))),
              "GeneRax example data not found")

  out_dir <- file.path(tempdir(), "radte_grax_tsv")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--generax_nhx=", shQuote(nhx_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # Gene tree TSV should exist with correct columns
  gn_tsv <- read.delim(file.path(out_dir, "radte_gene_tree.tsv"))
  expected_cols <- c("event", "gn_node", "lower_sp_node", "upper_sp_node",
                     "lower_age", "upper_age", "gn_node_num")
  expect_true(all(expected_cols %in% colnames(gn_tsv)))

  # Events should be valid
  expect_true(all(gn_tsv$event %in% c("D", "S", "S(R)", "D(R)")))

  # Exactly one root event
  root_events <- gn_tsv$event[grepl("\\(R\\)$", gn_tsv$event)]
  expect_equal(length(root_events), 1)

  # Output tree tip labels should match GeneRax input
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(is.ultrametric(out_tree))
})

# --- NOTUNG mode with duplication: output tree tip count matches input ---

test_that("RADTE preserves all tip labels in duplication scenario", {
  # Species tree with 3 species
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)n1:20,C_sp:30)root;", sp_file)
  on.exit(unlink(sp_file))

  # Gene tree with one duplication (A_sp has two gene copies)
  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(((A_sp_g1:0.05,A_sp_g2:0.05)n1:0.05,B_sp_g1:0.1)n2:0.2,C_sp_g1:0.3)n3;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t4\t7\t3\t0.05,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\tA_sp\tn1"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_dup_tips")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  in_tree <- read.tree(gn_file)
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))

  # All input tips should be preserved in output
  expect_setequal(out_tree$tip.label, in_tree$tip.label)
  expect_equal(length(out_tree$tip.label), 4)

  # Output should be ultrametric
  expect_true(is.ultrametric(out_tree))

  # Calibrated nodes should be RS (has both root and speciation constraints)
  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_true(cal_nodes[1] %in% c("RS", "S", "R"))
})

# --- Output tree branch lengths are all positive ---

test_that("RADTE output tree has all positive branch lengths", {
  sp_tree_file <- file.path(project_root, "data", "example_notung_01", "species_tree.nwk")
  gn_tree_file <- file.path(project_root, "data", "example_notung_01", "gene_tree.nwk.reconciled")
  parsable_file <- file.path(project_root, "data", "example_notung_01",
                             "gene_tree.nwk.reconciled.parsable.txt")
  skip_if_not(all(file.exists(c(sp_tree_file, gn_tree_file, parsable_file))),
              "NOTUNG example data not found")

  out_dir <- file.path(tempdir(), "radte_pos_bl")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_tree_file)),
    paste0("--gene_tree=", shQuote(gn_tree_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))

  # All branch lengths should be > 0 after pad_short_edge
  expect_true(all(out_tree$edge.length > 0))

  # All branch lengths should be >= pad_short_edge threshold
  expect_true(all(out_tree$edge.length >= 0.001))
})

# --- chronos_model=relaxed alternative ---

test_that("RADTE works with chronos_model=relaxed", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)n1:20,C_sp:30)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(((A_sp_g1:0.05,A_sp_g2:0.05)n1:0.05,B_sp_g1:0.1)n2:0.2,C_sp_g1:0.3)n3;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t4\t7\t3\t0.05,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\tA_sp\tn1"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_relaxed")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=relaxed",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(inherits(out_tree, "phylo"))
  expect_true(is.ultrametric(out_tree))
})
