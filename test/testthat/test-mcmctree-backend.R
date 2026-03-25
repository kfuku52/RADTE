radte_script <- file.path(project_root, "radte.r")

run_radte_mcmctree <- function(...) {
  out_dir <- tempfile("radte_mcmctree_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE))
  stderr_file <- tempfile()
  cmd <- paste("Rscript", shQuote(radte_script), paste(..., collapse = " "), "2>", shQuote(stderr_file))
  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)
  exit_code <- system(cmd, ignore.stdout = TRUE)
  stderr_output <- paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
  unlink(stderr_file)
  list(exit_code = exit_code, stderr = stderr_output)
}

test_that("MCMCTree calibration text widens fixed ages into hard bounds", {
  txt <- make_mcmctree_calibration_text(10, 10)
  expect_match(txt, "^B\\{")
  expect_false(grepl("B\\{10, 10\\}", txt))
})

test_that("MCMCTree root calibration text uses root-specific bound syntax", {
  txt <- make_mcmctree_root_calibration_text(20, 20)
  expect_match(txt, "^>")
  expect_match(txt, "<")
})

test_that("MCMCTree tree writer emits mirror labels for duplicated speciation nodes", {
  tree <- read.tree(text = "((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,(A_sp_g2:0.1,B_sp_g2:0.1)n2:0.1)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  n2_num <- get_node_num_by_name(tree, "n2")

  gn_node_table <- data.frame(
    event = c("S", "S", "D(R)"),
    gn_node = c("n1", "n2", "root"),
    gn_node_num = c(n1_num, n2_num, root_num),
    lower_sp_node = c("sp_ab", "sp_ab", "sp_root"),
    upper_sp_node = c("sp_ab", "sp_ab", NA),
    lower_age = c(5, 5, 20),
    upper_age = c(5, 5, 20),
    constraint_sp_node = c("sp_ab", "sp_ab", NA),
    shared_speciation_group = c("sp_ab", "sp_ab", NA),
    stringsAsFactors = FALSE
  )
  calibration_table <- data.frame(
    node = c(root_num, n1_num, n2_num),
    age.min = c(20, 5, 5),
    age.max = c(20, 5, 5),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  tree_info <- build_mcmctree_tree_text(tree, gn_node_table, calibration_table, root_num)
  expect_equal(tree_info$duplication_flag, 1L)
  expect_equal(length(gregexpr("#1", tree_info$tree_text, fixed = TRUE)[[1]]), 2)
  expect_equal(length(gregexpr("B{", tree_info$tree_text, fixed = TRUE)[[1]]), 1)
  expect_match(tree_info$tree_text, ">")
  expect_match(tree_info$tree_text, "<")
})

test_that("MCMCTree output parser extracts the posterior mean time tree", {
  out_file <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "Header",
      "Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels",
      "(A,B);",
      "(A:1.0,B:1.0);",
      "(A[&95%HPD={0.9,1.1}]:1.0,B:1.0);"
    ),
    out_file
  )
  on.exit(unlink(out_file))

  expect_equal(read_mcmctree_posterior_tree(out_file), "(A:1.0,B:1.0);")
})

test_that("MCMCTree control writer uses BDparas syntax accepted by current PAML", {
  ctl_file <- tempfile(fileext = ".ctl")
  on.exit(unlink(ctl_file))

  write_mcmctree_control_file(
    file = ctl_file,
    seqfile_name = "seqfile.phy",
    treefile_name = "input.trees"
  )

  ctl_lines <- readLines(ctl_file)
  expect_true(any(ctl_lines == "BDparas = 1 1 0.1 M"))
})

test_that("RADTE mcmctree backend requires seqfile", {
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

  result <- run_radte_mcmctree(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=100",
    "--dating_backend=mcmctree"
  )
  expect_true(result$exit_code != 0)
  expect_match(result$stderr, "mcmctree_seqfile", ignore.case = TRUE)
})
