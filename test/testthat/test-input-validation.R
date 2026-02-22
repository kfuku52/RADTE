# Input validation tests derived from GitHub issues #3, #4, #5

radte_script <- file.path(project_root, "radte.r")

# Helper: run RADTE with given args in a temp directory, return exit code and stderr
run_radte <- function(...) {
  out_dir <- file.path(tempdir(), paste0("radte_val_", as.integer(Sys.time())))
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
  list(exit_code = exit_code, stderr = stderr_output)
}

# --- Issue #3: Species tree with unlabeled internal nodes ---

test_that("RADTE errors on species tree with unlabeled internal nodes (issue #3)", {
  # Create a species tree with one unlabeled internal node (root has label, n1 does not)
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A:1,B:1):1,C:2)root;", sp_file)
  on.exit(unlink(sp_file))

  # Create a minimal gene tree and parsable file
  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "non-labeled", ignore.case = TRUE)
})

test_that("RADTE errors on gene tree with unlabeled internal nodes", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1):0.1,C_sp_g1:0.2);", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "gene tree.*non-labeled internal node", ignore.case = TRUE)
})

test_that("RADTE errors when all internal species tree nodes are unlabeled", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1):1,C_sp:2);", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "non-labeled", ignore.case = TRUE)
})

test_that("RADTE errors when species tree lacks branch lengths", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp,B_sp)n1,C_sp)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "species tree does not contain branch length", ignore.case = TRUE)
})

# --- Issue #4: Gene tree with polytomy ---

test_that("RADTE errors on gene tree with polytomy (issues #4, #5)", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  # Polytomy: 3 children at root
  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp_g1:0.1,B_sp_g1:0.1,C_sp_g1:0.1)n1;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "polytomy", ignore.case = TRUE)
})

test_that("RADTE errors on species tree with polytomy", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp:1,B_sp:1,C_sp:1)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "species tree.*polytomy", ignore.case = TRUE)
})

test_that("RADTE errors on species tree with singleton internal node", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1)n1:1,B_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp_g1:0.1,B_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "species tree.*polytomy", ignore.case = TRUE)
})

test_that("RADTE errors on gene tree with singleton internal node", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp:1,B_sp:1)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1)n1:0.1,B_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "input tree contains polytomy", ignore.case = TRUE)
})

test_that("RADTE errors when NOTUNG duplication annotation points to a tip node", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t4\t7\t3\t0.05,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tA_sp_g1\tA_sp\tn1"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "duplication event.*tip node", ignore.case = TRUE)
})

# --- Issue #4: Duplicate node names ---

test_that("check_gn_node_name_uniqueness detects duplicate node names (issue #4)", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,(C:1,D:1)n2:1)root;")

  # Unique labels are accepted
  gn_table <- data.frame(gn_node = c("n1", "n2"), stringsAsFactors = FALSE)
  expect_no_error(check_gn_node_name_uniqueness(gn_table, tree))

  # Duplicate internal labels should be rejected
  tree_dup <- tree
  tree_dup$node.label <- c("dup", "dup", "root")
  gn_table_dup <- data.frame(gn_node = c("dup"), stringsAsFactors = FALSE)
  expect_error(
    check_gn_node_name_uniqueness(gn_table_dup, tree_dup),
    "identical name|multiple nodes",
    ignore.case = TRUE
  )

  # Missing labels should also be rejected
  gn_table_bad <- data.frame(gn_node = c("missing_node"), stringsAsFactors = FALSE)
  expect_error(check_gn_node_name_uniqueness(gn_table_bad, tree),
               "does not contain node name|not found",
               ignore.case = TRUE)
})

test_that("RADTE errors on species tree with duplicated internal node labels", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,(C_sp:1,D_sp:1)n1:1)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)nA:0.1,(C_sp_g1:0.1,D_sp_g1:0.1)nB:0.1)nR;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "duplicated internal node label", ignore.case = TRUE)
})

# --- Issue #5: Zero branch lengths are padded ---

test_that("pad_branch_length handles zero branch lengths from gene trees (issue #5)", {
  tree <- read.tree(text = "((A:0,B:0.1)n1:0.05,C:0.2)root;")
  result <- pad_branch_length(tree, pad_size = 0.001)
  # The zero-length branch (A) should be padded
  a_idx <- which(result$edge[, 2] == which(result$tip.label == "A"))
  expect_gte(result$edge.length[a_idx], 0.001)
  # Non-zero branches should be unchanged
  c_idx <- which(result$edge[, 2] == which(result$tip.label == "C"))
  expect_equal(result$edge.length[c_idx], 0.2)
})

test_that("RADTE errors on gene trees with invalid tip label format", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((Aspg1:0.1,Bspg1:0.1)n1:0.1,Cspg1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "GENUS_SPECIES_GENEID", ignore.case = TRUE)
})

test_that("RADTE errors on gene trees with empty genus/species fields in tip labels", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A__g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "GENUS_SPECIES_GENEID", ignore.case = TRUE)
})

test_that("RADTE handles allS gene trees when GENEID contains underscores", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp:1,B_sp:1)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp_g_1:0.1,B_sp_g_1:0.1)n1;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_equal(result$exit_code, 0)
})

test_that("RADTE errors when gene tree contains species absent from species tree", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,X_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "were not found in the species tree", ignore.case = TRUE)
})

# --- Issue #4: Both mode flags specified ---

test_that("RADTE errors when both --notung_parsable and --generax_nhx are specified", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A:1,B:1)n1:1,C:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines("dummy", parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  nhx_file <- tempfile(fileext = ".nhx")
  writeLines("dummy", nhx_file)
  on.exit(unlink(nhx_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    paste0("--generax_nhx=", shQuote(nhx_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "Only one of", ignore.case = TRUE)
})

test_that("RADTE GeneRax mode errors when NHX lacks species annotation tag S", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  nhx_file <- tempfile(fileext = ".nhx")
  writeLines("((A_sp_g1:0.1[&&NHX:D=N],B_sp_g1:0.1[&&NHX:D=N])n1:0.1[&&NHX:D=N],C_sp_g1:0.2[&&NHX:D=N])n2:0.1[&&NHX:D=N];", nhx_file)
  on.exit(unlink(nhx_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--generax_nhx=", shQuote(nhx_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "required species annotation tag: S", ignore.case = TRUE)
})

test_that("RADTE GeneRax mode errors when NHX S annotation is absent from species tree", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  nhx_file <- tempfile(fileext = ".nhx")
  writeLines("((A_sp_g1:0.1[&&NHX:S=X_sp],B_sp_g1:0.1[&&NHX:S=B_sp])n1:0.1[&&NHX:S=n1],C_sp_g1:0.2[&&NHX:S=C_sp])n2:0.1[&&NHX:S=root];", nhx_file)
  on.exit(unlink(nhx_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--generax_nhx=", shQuote(nhx_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "lower_sp_node value\\(s\\) not found in the species tree", ignore.case = TRUE)
})

test_that("RADTE GeneRax mode works when NHX lacks duplication tag D", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  nhx_file <- tempfile(fileext = ".nhx")
  writeLines("((A_sp_g1:0.1[&&NHX:S=A_sp],B_sp_g1:0.1[&&NHX:S=B_sp])n1:0.1[&&NHX:S=n1],C_sp_g1:0.2[&&NHX:S=C_sp])n2:0.1[&&NHX:S=root];", nhx_file)
  on.exit(unlink(nhx_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--generax_nhx=", shQuote(nhx_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_equal(result$exit_code, 0)
})

test_that("RADTE GeneRax mode errors on unsupported duplication tag values", {
  skip_if_not(requireNamespace("treeio", quietly = TRUE), "treeio not installed")

  sp_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp:1,B_sp:1)root;", sp_file)
  on.exit(unlink(sp_file))

  nhx_file <- tempfile(fileext = ".nhx")
  writeLines(
    "((A_sp_g1:0.1[&&NHX:S=A_sp:D=maybe],A_sp_g2:0.1[&&NHX:S=A_sp:D=maybe])n1:0.1[&&NHX:S=A_sp:D=maybe],B_sp_g1:0.2[&&NHX:S=B_sp:D=maybe])n2:0.1[&&NHX:S=root:D=maybe];",
    nhx_file
  )
  on.exit(unlink(nhx_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--generax_nhx=", shQuote(nhx_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "unsupported duplication tag value", ignore.case = TRUE)
})

test_that("RADTE errors when duplicated-mode gene tree lacks branch lengths", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)n1:20,C_sp:30)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(((A_sp_g1,A_sp_g2)n1,B_sp_g1)n2,C_sp_g1)n3;", gn_file)
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

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "gene tree for chronos does not contain branch length", ignore.case = TRUE)
})

test_that("RADTE NOTUNG mode errors clearly when lower/upper age bounds remain unresolved", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn2\tX_sp\troot"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "unresolved age bound", ignore.case = TRUE)
})

# --- Neither mode flag specified ---

test_that("RADTE errors when neither --notung_parsable nor --generax_nhx is specified", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A:1,B:1)n1:1,C:2)root;", sp_file)
  on.exit(unlink(sp_file))

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "should be specified", ignore.case = TRUE)
})

test_that("RADTE with a single argument reports missing mode flags", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("(A_sp:1,B_sp:1)root;", sp_file)
  on.exit(unlink(sp_file))

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file))
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "should be specified", ignore.case = TRUE)
  expect_false(grepl("cannot change working directory", result$stderr, ignore.case = TRUE))
})

test_that("RADTE errors with explicit message when --species_tree is missing", {
  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file))

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "Missing required argument\\(s\\): --species_tree", ignore.case = TRUE)
})

test_that("RADTE errors with explicit message when --gene_tree is missing in NOTUNG mode", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "--gene_tree should be specified in Notung mode", ignore.case = TRUE)
})

test_that("RADTE errors on non-numeric --max_age", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=abc", "--chronos_lambda=1", "--chronos_model=discrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "--max_age should be a positive finite number", ignore.case = TRUE)
})

test_that("RADTE errors on non-numeric --chronos_lambda", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=abc", "--chronos_model=discrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "--chronos_lambda should be a non-negative finite number", ignore.case = TRUE)
})

test_that("RADTE errors when NOTUNG parsable has duplicated gn_node rows", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn2\tn1\troot",
    "#D\tn2\tn1\troot"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "duplicated gn_node", ignore.case = TRUE)
})

test_that("RADTE errors on invalid --allow_constraint_drop", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--allow_constraint_drop=maybe"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "--allow_constraint_drop should be boolean", ignore.case = TRUE)
})

test_that("RADTE errors on unsupported --chronos_model values", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=invalid_model"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "--chronos_model should be one of", ignore.case = TRUE)
})

test_that("RADTE suggests the correct model name for --chronos_model=difscrete", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=difscrete"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "Did you mean \"discrete\"", fixed = TRUE)
})

test_that("RADTE errors on invalid --chronos_attempt_timeout_sec", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--chronos_attempt_timeout_sec=-1"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "--chronos_attempt_timeout_sec should be a non-negative number", ignore.case = TRUE)
})

test_that("RADTE errors on invalid --chronos_total_timeout_sec", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:1,B_sp:1)n1:1,C_sp:2)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,C_sp_g1:0.2)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100", "--chronos_lambda=1", "--chronos_model=discrete",
    "--chronos_total_timeout_sec=not_a_number"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "--chronos_total_timeout_sec should be a non-negative number", ignore.case = TRUE)
})
