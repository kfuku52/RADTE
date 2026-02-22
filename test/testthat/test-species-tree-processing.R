# Tests for species tree processing logic (quoted numeric labels, PLACEHOLDER substitution)

radte_script <- file.path(project_root, "radte.r")

# Helper: run RADTE and return exit code + stderr
run_radte_sp <- function(...) {
  out_dir <- file.path(tempdir(), paste0("radte_sp_", as.integer(Sys.time())))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))
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

# --- Species tree with quoted numeric node labels ---

test_that("RADTE handles species tree with quoted numeric node labels (line 376)", {
  # Species trees from some tools have numeric labels quoted like '123'
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((Homo_sapiens:10,Mus_musculus:10)'100':20,Danio_rerio:30)'200';", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((Homo_sapiens_g1:0.1,Mus_musculus_g1:0.1)n1:0.2,Danio_rerio_g1:0.3)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t3\t5\t3\t0.1,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_quoted_labels")
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

  # Verify species tree TSV has unquoted numeric labels
  sp_tsv <- read.delim(file.path(out_dir, "radte_species_tree.tsv"))
  # The node labels '100' and '200' should have quotes stripped
  expect_true("100" %in% sp_tsv$node | "200" %in% sp_tsv$node)
})

test_that("RADTE preserves non-temporary PLACEHOLDER substring in species node labels", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)XPLACEHOLDER100:20,C_sp:30)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.2,C_sp_g1:0.3)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t3\t5\t3\t0.1,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_placeholder_preserve")
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

  sp_tsv <- read.delim(file.path(out_dir, "radte_species_tree.tsv"))
  expect_true("XPLACEHOLDER100" %in% sp_tsv$node)
})

# --- Species tree with all species in gene tree (no drop needed in allS) ---

test_that("RADTE allS mode works when gene tree covers all species", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)n1:20,C_sp:30)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.2,C_sp_g1:0.3)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t3\t5\t3\t0.1,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_allS_full")
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

  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_equal(cal_nodes[1], "allS")

  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_equal(length(out_tree$tip.label), 3)
  expect_true(is.ultrametric(out_tree))
})

# --- allS mode with species subset (gene tree has fewer species) ---

test_that("RADTE allS mode drops extra species from species tree", {
  # 4 species in species tree, but only 3 in gene tree (D_sp missing)
  # Need >=3 tips so there's a non-root speciation node for calibration_table_S
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)n1:10,(C_sp:10,D_sp:10)n2:10)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.2,C_sp_g1:0.3)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t3\t5\t4\t0.1,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_allS_subset")
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

  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_equal(cal_nodes[1], "allS")

  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  # Output should only have 3 tips (the gene tree species, D_sp dropped)
  expect_equal(length(out_tree$tip.label), 3)
})

test_that("RADTE allS mode works for two-species trees (root-only speciation)", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp:10,B_sp:10)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp_g1:0.1,B_sp_g1:0.1)n1;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t2\t3\t2\t0.1,0.1\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), "radte_allS_two_tip")
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

  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_equal(cal_nodes[1], "allS")

  cal_used <- read.delim(file.path(out_dir, "radte_calibration_used.tsv"))
  expect_equal(nrow(cal_used), 1)
  expect_equal(cal_used$event[1], "R")
  expect_equal(cal_used$age.min[1], 10)
  expect_equal(cal_used$age.max[1], 10)

  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_equal(length(out_tree$tip.label), 2)
  expect_true(is.ultrametric(out_tree))
})
