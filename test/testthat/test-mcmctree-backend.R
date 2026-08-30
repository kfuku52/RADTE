radte_script <- file.path(project_root, "radte.r")

run_radte_mcmctree <- function(...) {
  out_dir <- tempfile("radte_mcmctree_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE))
  stderr_file <- tempfile()
  file.create(stderr_file)
  args <- c(...)
  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)
  exit_code <- system2(
    "Rscript",
    c(shQuote(radte_script), args),
    stdout = FALSE,
    stderr = stderr_file
  )
  stderr_output <- paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
  unlink(stderr_file)
  list(exit_code = exit_code, stderr = stderr_output)
}

test_that("MCMCTree calibration text widens fixed ages into narrow PAML prior intervals", {
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
  expect_equal(sort(tree_info$mirror_table$mirror_role), c("driver", "mirror"))
  expect_equal(sum(tree_info$mirror_table$prior_emitted), 1)
})

test_that("MCMCTree calibration validation preserves feasible overlapping bounds", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  calibration_table <- data.frame(
    node = c(root_num, n1_num),
    age.min = c(5, 6),
    age.max = c(12, 9),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  original <- calibration_table

  expect_silent(validate_mcmctree_calibration_constraints(tree, calibration_table))
  expect_identical(calibration_table, original)
})

test_that("MCMCTree calibration validation rejects impossible temporal order", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  calibration_table <- data.frame(
    node = c(root_num, n1_num),
    age.min = c(5, 12),
    age.max = c(12, 14),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  expect_error(
    validate_mcmctree_calibration_constraints(tree, calibration_table),
    "temporally infeasible"
  )
})

test_that("MCMCTree mirror writer rejects missing member calibration", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,(C:1,D:1)n2:1)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  n2_num <- get_node_num_by_name(tree, "n2")
  gn_node_table <- data.frame(
    event = c("S", "S", "S(R)"),
    gn_node = c("n1", "n2", "root"),
    gn_node_num = c(n1_num, n2_num, root_num),
    shared_speciation_group = c("sp_ab", "sp_ab", NA),
    stringsAsFactors = FALSE
  )
  calibration_table <- data.frame(
    node = c(root_num, n1_num),
    age.min = c(10, 5),
    age.max = c(12, 7),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  expect_error(
    build_mcmctree_tree_text(tree, gn_node_table, calibration_table, root_num),
    "missing calibration row"
  )
})

test_that("MCMCTree mirror writer rejects inconsistent member bounds", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,(C:1,D:1)n2:1)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  n2_num <- get_node_num_by_name(tree, "n2")
  gn_node_table <- data.frame(
    event = c("S", "S", "S(R)"),
    gn_node = c("n1", "n2", "root"),
    gn_node_num = c(n1_num, n2_num, root_num),
    shared_speciation_group = c("sp_ab", "sp_ab", NA),
    stringsAsFactors = FALSE
  )
  calibration_table <- data.frame(
    node = c(root_num, n1_num, n2_num),
    age.min = c(10, 5, 6),
    age.max = c(12, 7, 7),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  expect_error(
    build_mcmctree_tree_text(tree, gn_node_table, calibration_table, root_num),
    "inconsistent calibration bounds"
  )
})

test_that("shared speciation posterior QA reports equal ages", {
  tree <- read.tree(text = "((A:5,B:5)n1:5,(C:5,D:5)n2:5)root;")
  gn_node_table <- data.frame(
    event = c("S", "S", "S(R)"),
    gn_node = c("n1", "n2", "root"),
    gn_node_num = c(
      get_node_num_by_name(tree, "n1"),
      get_node_num_by_name(tree, "n2"),
      get_root_num(tree)
    ),
    shared_speciation_group = c("sp_ab", "sp_ab", NA),
    stringsAsFactors = FALSE
  )

  summary_table <- summarize_shared_speciation_ages(tree, gn_node_table)
  expect_equal(nrow(summary_table), 1)
  expect_equal(summary_table$species_node, "sp_ab")
  expect_equal(summary_table$member_count, 2)
  expect_equal(summary_table$max_member_age_diff, 0)
})

test_that("shared speciation posterior QA rejects unequal ages", {
  tree <- read.tree(text = "((A:5,B:5)n1:5,(C:4,D:4)n2:6)root;")
  gn_node_table <- data.frame(
    event = c("S", "S", "S(R)"),
    gn_node = c("n1", "n2", "root"),
    gn_node_num = c(
      get_node_num_by_name(tree, "n1"),
      get_node_num_by_name(tree, "n2"),
      get_root_num(tree)
    ),
    shared_speciation_group = c("sp_ab", "sp_ab", NA),
    stringsAsFactors = FALSE
  )

  expect_error(
    summarize_shared_speciation_ages(tree, gn_node_table, tolerance = 1e-8),
    "posterior ages differ"
  )
})

test_that("calibration output identifies mirror roles and emitted prior", {
  calibration_table <- data.frame(
    node = c(7L, 8L, 9L, 10L),
    age.min = c(10, 5, 5, 2),
    age.max = c(12, 7, 7, 4),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  gn_node_table <- data.frame(
    gn_node_num = c(7L, 8L, 9L, 10L),
    gn_node = c("root", "n1", "n2", "d1"),
    event = c("S(R)", "S", "S", "D"),
    shared_speciation_group = c(NA, "sp_ab", "sp_ab", NA),
    stringsAsFactors = FALSE
  )
  mirror_table <- data.frame(
    shared_speciation_group = c("sp_ab", "sp_ab"),
    node = c(8L, 9L),
    gn_node = c("n1", "n2"),
    mirror_label = c("#1", "#1"),
    mirror_role = c("driver", "mirror"),
    prior_emitted = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  output <- annotate_calibration_output(
    calibration_table,
    gn_node_table,
    mirror_table = mirror_table,
    emitted_nodes = c(7L, 8L, 9L)
  )
  expect_equal(output$mirror_role[output$node == 8L], "driver")
  expect_equal(output$mirror_role[output$node == 9L], "mirror")
  expect_true(output$prior_emitted[output$node == 8L])
  expect_false(output$prior_emitted[output$node == 9L])
  expect_false(output$prior_emitted[output$node == 10L])
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
    treefile_name = "input.trees",
    seed = 41L
  )

  ctl_lines <- readLines(ctl_file)
  expect_true(any(ctl_lines == "BDparas = 1 1 0.1 M"))
  expect_true(any(ctl_lines == "seed = 41"))
})

test_that("MCMCTree external execution honors its wall-time limit", {
  skip_on_os("windows")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  root_num <- get_root_num(tree)
  calibration_table <- data.frame(
    node = root_num,
    age.min = 2,
    age.max = 3,
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  seqfile <- tempfile(fileext = ".phy")
  fake_bin <- tempfile(pattern = "fake_mcmctree_")
  workdir <- tempfile(pattern = "mcmctree_timeout_")
  writeLines("3 1", seqfile)
  writeLines(c("#!/usr/bin/env Rscript", "Sys.sleep(5)"), fake_bin)
  Sys.chmod(fake_bin, "0755")
  on.exit(unlink(c(seqfile, fake_bin, workdir), recursive = TRUE))

  expect_error(
    run_mcmctree_backend(
      phy = tree,
      gn_node_table = data.frame(),
      calibration_table = calibration_table,
      root_num = root_num,
      seqfile = seqfile,
      bin = fake_bin,
      workdir = workdir,
      timeout_sec = 0.1
    ),
    "timeout|exit code 124",
    ignore.case = TRUE
  )
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

test_that("PAML integration preserves bounds and verifies mirror ages", {
  run_external <- tolower(Sys.getenv("RADTE_RUN_PAML_TESTS", unset = "false")) %in%
    c("1", "true", "yes")
  skip_if_not(run_external, "Set RADTE_RUN_PAML_TESTS=true to run the external PAML test")
  mcmctree_bin <- Sys.which("mcmctree")
  if (!nzchar(mcmctree_bin)) {
    stop("RADTE_RUN_PAML_TESTS=true but mcmctree is not available in PATH")
  }

  fixture_dir <- tempfile("radte_paml_fixture_")
  out_dir <- tempfile("radte_paml_output_")
  dir.create(fixture_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(c(fixture_dir, out_dir), recursive = TRUE))

  sp_file <- file.path(fixture_dir, "species.nwk")
  gn_file <- file.path(fixture_dir, "gene.nwk")
  parsable_file <- file.path(fixture_dir, "notung.txt")
  bounds_file <- file.path(fixture_dir, "bounds.tsv")
  seq_file <- file.path(fixture_dir, "alignment.phy")

  writeLines("((A_sp:8,B_sp:8)sp_ab:2,C_sp:10)sp_root;", sp_file)
  writeLines(
    "(((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,(A_sp_g2:0.1,B_sp_g2:0.1)n2:0.1)d1:0.1,C_sp_g1:0.3)groot;",
    gn_file
  )
  writeLines(
    c(
      "1\t0\t0\t0\t5\t9\t3\t0.1,1.0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
      "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
      "",
      "#D\tDuplication\tL.Bound\tU.Bound",
      "#D\td1\tsp_ab\tsp_root"
    ),
    parsable_file
  )
  writeLines(
    c(
      "node\tage_min\tage_max",
      "sp_ab\t6\t9",
      "sp_root\t5\t12"
    ),
    bounds_file
  )
  sequences <- c(
    A_sp_g1 = paste(rep("ACGT", 30), collapse = ""),
    B_sp_g1 = paste(rep("ACGA", 30), collapse = ""),
    A_sp_g2 = paste(rep("AGGT", 30), collapse = ""),
    B_sp_g2 = paste(rep("AGGA", 30), collapse = ""),
    C_sp_g1 = paste(rep("TCGT", 30), collapse = "")
  )
  writeLines(
    c(
      paste(length(sequences), nchar(sequences[[1]])),
      paste(names(sequences), sequences, sep = "  ")
    ),
    seq_file
  )

  rscript_bin <- file.path(R.home("bin"), "Rscript")
  args <- c(
    radte_script,
    paste0("--species_tree=", sp_file),
    paste0("--gene_tree=", gn_file),
    paste0("--notung_parsable=", parsable_file),
    paste0("--species_node_bounds_tsv=", bounds_file),
    "--max_age=20",
    "--dating_backend=mcmctree",
    paste0("--mcmctree_seqfile=", seq_file),
    paste0("--mcmctree_bin=", mcmctree_bin),
    paste0("--mcmctree_workdir=", file.path(out_dir, "work")),
    "--mcmctree_burnin=20",
    "--mcmctree_sampfreq=1",
    "--mcmctree_nsample=50"
  )
  stdout_file <- file.path(out_dir, "radte.stdout.log")
  stderr_file <- file.path(out_dir, "radte.stderr.log")
  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)
  exit_code <- system2(
    rscript_bin,
    args = args,
    stdout = stdout_file,
    stderr = stderr_file
  )
  if (exit_code != 0) {
    diagnostics <- paste(
      c(readLines(stdout_file, warn = FALSE), readLines(stderr_file, warn = FALSE)),
      collapse = "\n"
    )
    fail(paste("RADTE/PAML integration failed:\n", diagnostics))
    return(invisible(NULL))
  }

  gene_table <- read.delim(file.path(out_dir, "radte_gene_tree.tsv"))
  mirrored <- gene_table[gene_table$gn_node %in% c("n1", "n2"), ]
  expect_equal(sort(mirrored$lower_age), c(6, 6))
  expect_equal(sort(mirrored$upper_age), c(9, 9))
  expect_true(all(mirrored$shared_speciation_group == "sp_ab"))

  calibration_used <- read.delim(file.path(out_dir, "radte_calibration_used.tsv"))
  mirror_calibration <- calibration_used[
    calibration_used$shared_speciation_group == "sp_ab" &
      !is.na(calibration_used$shared_speciation_group),
  ]
  expect_equal(sort(mirror_calibration$mirror_role), c("driver", "mirror"))
  expect_equal(sum(mirror_calibration$prior_emitted), 1)

  ages <- read.delim(file.path(out_dir, "radte_shared_speciation_ages.tsv"))
  expect_equal(ages$species_node, "sp_ab")
  expect_equal(ages$member_count, 2)
  expect_lte(ages$max_member_age_diff, ages$tolerance)

  input_tree_text <- paste(
    readLines(file.path(out_dir, "radte_mcmctree_input.trees"), warn = FALSE),
    collapse = ""
  )
  expect_equal(length(gregexpr("#1", input_tree_text, fixed = TRUE)[[1]]), 2)
  expect_equal(length(gregexpr("B{", input_tree_text, fixed = TRUE)[[1]]), 1)
})
