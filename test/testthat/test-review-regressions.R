test_that("map and bounds TSVs preserve numeric-looking IDs and uppercase headers", {
  map <- tempfile()
  bounds <- tempfile()
  on.exit(unlink(c(map, bounds)))
  writeLines(c("LABEL\tSPECIES", "001\t101", "002\t102"), map)
  writeLines(c("NODE\tAGE_MIN\tAGE_MAX", "001\t8\t10"), bounds)
  expect_identical(read_species_map_tsv(map)$label, c("001", "002"))
  expect_identical(read_species_node_bounds_tsv(bounds)$node, "001")
  expect_equal(read_species_node_bounds_tsv(bounds)$age_min, 8)
})

test_that("feasible overlapping age ranges remain unchanged and constrain the result", {
  tree <- read.tree(text = "((A:9,B:9)child:1,C:10)root;")
  nodes <- data.frame(gn_node_num = c(4L, 5L), lower_age = c(5, 8), upper_age = c(15, 10))
  expect_identical(stabilize_descendant_constraints(nodes, tree, 4L)$gn_node_table, nodes)
  expect_equal(nrow(find_descendant_constraint_conflicts(nodes, tree, 4L)), 0)
  calibration <- data.frame(node = c(4L, 5L), age.min = c(5, 8), age.max = c(15, 10))
  dated <- build_dated_tree_without_chronos(tree, calibration, 4L)
  expect_silent(validate_dated_tree(dated, tree, calibration))
  original <- calibration
  calibration$age.min <- calibration$age.max <- c(10, 15)
  expect_s3_class(build_dated_tree_without_chronos(tree, calibration, 4L), "try-error")
  expect_equal(original$age.min, c(5, 8))
  wrong <- tree
  wrong$edge.length <- c(5, 5, 5, 10)
  expect_error(validate_dated_tree(wrong, tree, original), "original calibration bounds")
})

test_that("finite-timeout fork and inline execution inherit the same random stream", {
  skip_on_os("windows")
  runner <- run_chronos_once
  environment(runner) <- list2env(list(radte_chronos = function(...) runif(6)),
                                 parent = environment(runner))
  invoke <- function(timeout) {
    set.seed(42)
    runner(NULL, NULL, NULL, 1, "discrete", timeout)
  }
  expect_identical(invoke(5), invoke(5))
  expect_identical(invoke(5), invoke(Inf))
})

test_that("real chronos repeats with finite timeout and restores the caller RNG", {
  tree <- read.tree(text = "((A:0.1,B:0.2)x:0.3,(C:0.5,D:0.4)y:0.2)root;")
  calibration <- data.frame(node = 5L, age.min = 10, age.max = 20, soft.bounds = FALSE)
  set.seed(71)
  saved <- .Random.seed
  run <- function() suppressWarnings(run_chronos_with_restarts(
    tree, calibration, chronos.control(), 1, "discrete", max_restarts = 1,
    seed_base = 42, attempt_timeout_sec = 5
  ))
  one <- run()
  two <- run()
  expect_true(one$success)
  expect_true(two$success)
  expect_identical(one$chronos_out$edge.length, two$chronos_out$edge.length)
  expect_identical(.Random.seed, saved)
})

test_that("allS transfers ages without changing gene node numbers or topology", {
  species <- read.tree(text = "((C_sp:7,D_sp:7)sp_cd:3,(A_sp:5,B_sp:5)sp_ab:5)sp_root;")
  gene <- read.tree(text = "((A_sp_g1:1,B_sp_g1:1)gn_ab:1,(C_sp_g1:1,D_sp_g1:1)gn_cd:1)gn_root;")
  out <- transfer_species_ages(gene, species, build_species_parser("legacy"))
  expect_identical(out$edge, gene$edge)
  expect_identical(out$tip.label, gene$tip.label)
  expect_identical(out$node.label, gene$node.label)
  ages <- max(node.depth.edgelength(out)) - node.depth.edgelength(out)
  expect_equal(ages[c(5, 6, 7)], c(10, 5, 7))
  discordant <- read.tree(text = "((A_sp_g1:1,C_sp_g1:1)ac:1,(B_sp_g1:1,D_sp_g1:1)bd:1)root;")
  expect_error(transfer_species_ages(discordant, species, build_species_parser("legacy")), "topology")
})

test_that("GeneRax and NOTUNG resolve duplication bounds through the same species nodes", {
  skip_if_not_installed("treeio")
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  writeLines("((A_sp:5,B_sp:5)sp_ab:5,C_sp:10)sp_root;", file.path(dir, "sp.nwk"))
  writeLines("((A_sp_g1:1,A_sp_g2:1)dup:1,B_sp_g1:2)gr;", file.path(dir, "gene.nwk"))
  writeLines("#D dup A_sp sp_ab", file.path(dir, "events"))
  writeLines(paste0("((A_sp_g1:1[&&NHX:S=A_sp:D=N],A_sp_g2:1[&&NHX:S=A_sp:D=N])",
                   "dup:1[&&NHX:S=A_sp:D=Y],B_sp_g1:2[&&NHX:S=B_sp:D=N])gr[&&NHX:S=sp_ab:D=N];"),
             file.path(dir, "gene.nhx"))
  args <- list(species_tree = file.path(dir, "sp.nwk"), gene_tree = file.path(dir, "gene.nwk"),
               notung_parsable = file.path(dir, "events"), generax_nhx = file.path(dir, "gene.nhx"), max_age = 20)
  parser <- build_species_parser("legacy")
  sp <- read_radte_species_context(args, parser)
  notung <- read_radte_gene_context(args, "notung", sp$tree, sp$nodes, parser)$nodes
  generax <- read_radte_gene_context(args, "generax", sp$tree, sp$nodes, parser)$nodes
  cols <- c("lower_sp_node", "upper_sp_node", "lower_age", "upper_age")
  one <- notung[notung$gn_node == "dup", cols]
  two <- generax[generax$gn_node == "dup", cols]
  rownames(one) <- rownames(two) <- NULL
  expect_identical(one, two)
  expect_equal(c(one$lower_age, one$upper_age), c(0, 5))
  writeLines("#D dup A_sp TYPO_UNKNOWN", args$notung_parsable)
  expect_error(read_radte_gene_context(args, "notung", sp$tree, sp$nodes, parser), "not found.*TYPO_UNKNOWN")
  writeLines("#D dup A_sp NA", args$notung_parsable)
  expect_error(read_radte_gene_context(args, "notung", sp$tree, sp$nodes, parser), "unresolved age")
})

test_that("staging refuses aliases of the original alignment without changing its hash", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  input <- file.path(dir, "seqfile.phy")
  writeLines(c("2 1", "A A", "B A"), input)
  before <- compute_file_sha256(input)
  expect_error(stage_external_input_file(input, input), "collision")
  expect_error(stage_external_input_file(input, file.path(dir, ".", "seqfile.phy")), "collision")
  if (.Platform$OS.type == "unix") {
    alias <- file.path(dir, "alias")
    expect_true(file.symlink(input, alias))
    expect_error(stage_external_input_file(alias, input), "collision")
    expect_error(stage_external_input_file(input, alias), "collision")
  }
  expect_identical(compute_file_sha256(input), before)
  snapshot <- file.path(dir, "copy.phy")
  stage_external_input_file(input, snapshot)
  writeLines("modified snapshot", snapshot)
  expect_identical(compute_file_sha256(input), before)
})

test_that("MCMCTree checks collisions and calibration before cleaning its work directory", {
  skip_on_os("windows")
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  tree <- read.tree(text = "(A:1,B:1)root;")
  input <- file.path(dir, "seqfile.phy")
  old_output <- file.path(dir, "out.txt")
  writeLines("alignment", input)
  writeLines("previous result", old_output)
  run <- function(cal) run_mcmctree_backend(tree, data.frame(), cal, 3L, input, bin = Sys.which("true"), workdir = dir)
  cal <- data.frame(node = 3L, age.min = 1, age.max = 2)
  expect_error(run(cal), "collision")
  expect_identical(readLines(input), "alignment")
  expect_identical(readLines(old_output), "previous result")
  cal$age.max <- 0
  expect_error(run(cal), "younger than")
  expect_identical(readLines(old_output), "previous result")
})

test_that("MCMCTree supports relative executable paths and rejects malformed posterior trees", {
  skip_on_os("windows")
  dir <- tempfile()
  dir.create(dir)
  old <- getwd()
  on.exit({
    setwd(old)
    unlink(dir, recursive = TRUE)
  })
  setwd(dir)
  dir.create("relative bin")
  input <- file.path(dir, "alignment")
  writeLines("2 1", input)
  bin <- file.path("relative bin", "sentinel")
  writeLines(c("#!/bin/sh", "exit 7"), bin)
  Sys.chmod(bin, "0755")
  tree <- read.tree(text = "(A:1,B:1)root;")
  cal <- data.frame(node = 3L, age.min = 1, age.max = 2)
  run <- function(binary = bin) run_mcmctree_backend(tree, data.frame(), cal, 3L, input, bin = binary, workdir = "work")
  expect_error(run(), "exit code 7")
  expect_error(run(normalizePath(bin)), "exit code 7")
  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(normalizePath("relative bin"), old_path, sep = .Platform$path.sep))
  expect_error(run("sentinel"), "exit code 7")
  for (posterior in c("(A:-1,B:2);", "(A:1,B:2);", "(A:1,C:1);")) {
    writeLines(c("#!/bin/sh", "cat > out.txt <<'TREE'", "Species tree for FigTree", "(A:1,B:1);", posterior, "TREE"), bin)
    expect_error(run(), "Invalid dated tree")
  }
  writeLines(c("#!/bin/sh", "cat > out.txt <<'TREE'", "Species tree for FigTree", "(B:1,A:1);", "(B:1,A:1);", "TREE"), bin)
  result <- run()
  expect_identical(result$tree$edge, tree$edge)
  expect_identical(result$tree$tip.label, tree$tip.label)
})

test_that("MCMCTree writer handles deep trees without recursion", {
  tree <- stree(4000, type = "left")
  tree$node.label <- paste0("n", seq_len(tree$Nnode))
  root <- get_root_num(tree)
  cal <- data.frame(node = root, age.min = 1, age.max = 2)
  out <- build_mcmctree_tree_text(tree, data.frame(), cal, root)
  expect_match(out$tree_text, ";$")
  parsed <- read.tree(text = out$tree_text)
  expect_identical(parsed$tip.label, tree$tip.label)
  expect_false(anyNA(match_internal_nodes_by_clade(tree, parsed)))
})

test_that("dependency installer accepts final-attempt success and rejects total failure", {
  expressions <- parse(file.path(project_root, "test", "install_dependencies.R"))
  definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                          identical(x[[2]], as.name("install_cran_with_retry")), as.list(expressions))[[1]]
  for (success_at in c(0L, 1L, 3L, 4L)) {
    calls <- 0L
    env <- new.env(parent = baseenv())
    env$requireNamespace <- function(...) calls >= success_at
    env$install.packages <- function(...) { calls <<- calls + 1L }
    eval(definition, env)
    if (success_at <= 3L) {
      expect_true(env$install_cran_with_retry("stub"))
      expect_equal(calls, success_at)
    } else {
      expect_error(env$install_cran_with_retry("stub"), "after 3 attempts")
      expect_equal(calls, 3L)
    }
  }
})

test_that("a late PDF failure leaves the previous complete output generation intact", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  sp <- file.path(dir, "sp.nwk")
  gene <- file.path(dir, "gene.nwk")
  events <- file.path(dir, "events")
  writeLines("(A_sp:10,B_sp:10)sp_root;", sp)
  writeLines("(A_sp_g1:1,B_sp_g1:1)gn_root;", gene)
  writeLines(character(), events)
  args <- c(paste0("--species_tree=", sp), paste0("--gene_tree=", gene),
            paste0("--notung_parsable=", events), paste0("--outdir=", dir),
            "--max_age=20", "--chronos_lambda=1", "--chronos_model=discrete")
  result <- radte_main(args)
  outputs <- list.files(dir, pattern = "^radte_", full.names = TRUE)
  hashes <- tools::md5sum(outputs)
  writer <- write_radte_outputs
  environment(writer) <- list2env(list(atomic_save_tree_pdf = function(...) stop("injected PDF failure")),
                                  parent = environment(writer))
  runner <- radte_main
  environment(runner) <- list2env(list(write_radte_outputs = writer), parent = environment(runner))
  expect_error(runner(args), "injected PDF failure")
  expect_identical(tools::md5sum(outputs), hashes)
  expect_false(dir.exists(file.path(dir, ".radte.radte-lock")))
  manifest <- read.delim(result$manifest, colClasses = "character")
  expect_equal(manifest$value[manifest$key == "run_status"], "complete")
  expect_equal(manifest$value[manifest$key == "option.seed"], "1")
  expect_true("rng_kind" %in% manifest$key)
  expect_true(any(startsWith(manifest$key, "output_sha256.")))
})

test_that("publication failure rolls back all artifacts and excludes concurrent writers", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  tree <- file.path(dir, "x_gene_tree_output.nwk")
  manifest <- file.path(dir, "x_run_manifest.tsv")
  writeLines("old tree", tree)
  writeLines("old manifest", manifest)
  transaction <- begin_output_transaction(dir, "x")
  on.exit(cleanup_output_transaction(transaction), add = TRUE)
  expect_error(begin_output_transaction(dir, "x"), "Another run owns")
  writeLines("new tree", file.path(transaction$stage, basename(tree)))
  writeLines("new manifest", file.path(transaction$stage, basename(manifest)))
  publish <- commit_output_transaction
  environment(publish) <- list2env(list(file.rename = function(from, to) {
    if (from == file.path(transaction$stage, basename(manifest))) return(FALSE)
    base::file.rename(from, to)
  }), parent = environment(publish))
  expect_error(publish(transaction, basename(manifest)), "Could not publish")
  expect_identical(readLines(tree), "old tree")
  expect_identical(readLines(manifest), "old manifest")
})

test_that("minimum-edge projection preserves feasible ages and rejects hard-bound conflicts", {
  tree <- read.tree(text = "((A:1,B:1)child:0.0001,C:1.0001)root;")
  cal <- data.frame(node = 4L, age.min = 1, age.max = 2, soft.bounds = FALSE)
  out <- build_dated_tree_without_chronos(tree, cal, 4L, min_edge = 0.001, require_root = FALSE)
  expect_s3_class(out, "phylo")
  expect_equal(max(node.depth.edgelength(out)), 1.0001)
  expect_true(all(out$edge.length >= 0.001 - 1e-12))
  cal <- rbind(cal, data.frame(node = 5L, age.min = 1, age.max = 1, soft.bounds = FALSE))
  cal$age.min[[1]] <- cal$age.max[[1]] <- 1.0001
  expect_s3_class(build_dated_tree_without_chronos(tree, cal, 4L, min_edge = 0.001), "try-error")
})
