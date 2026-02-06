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

# --- Issue #4: Duplicate node names ---

test_that("check_gn_node_name_uniqueness detects duplicate node names (issue #4)", {
  # Tree with two nodes named "dup"
  # We cannot easily create this with read.tree, so test the function directly
  tree <- read.tree(text = "((A:1,B:1)n1:1,(C:1,D:1)n2:1)root;")
  # Simulate a table referencing a valid node
  gn_table <- data.frame(gn_node = c("n1", "n2"), stringsAsFactors = FALSE)
  expect_no_error(check_gn_node_name_uniqueness(gn_table, tree))

  # Simulate a table referencing a node that doesn't exist (as if duplicated and removed)
  gn_table_bad <- data.frame(gn_node = c("missing_node"), stringsAsFactors = FALSE)
  expect_error(check_gn_node_name_uniqueness(gn_table_bad, tree),
               "multiple nodes|not found|identical name",
               ignore.case = TRUE)
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
