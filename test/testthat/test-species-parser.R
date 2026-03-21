radte_script <- file.path(project_root, "radte.r")

run_radte_species_parser <- function(...) {
  out_dir <- tempfile(pattern = "radte_species_parser_")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE))
  stderr_file <- tempfile()
  args <- paste(...)
  cmd <- paste("Rscript", shQuote(radte_script), args, "2>", shQuote(stderr_file))
  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)
  exit_code <- system(cmd, ignore.stdout = TRUE)
  stderr_output <- paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
  stdout_files <- list.files(out_dir, full.names = TRUE)
  output_tree <- NULL
  output_tree_path <- file.path(out_dir, "radte_gene_tree_output.nwk")
  if (file.exists(output_tree_path)) {
    output_tree <- read.tree(output_tree_path)
  }
  unlink(stderr_file)
  list(exit_code = exit_code, stderr = stderr_output, out_dir = out_dir, output_tree = output_tree)
}

test_that("species parser defaults to legacy behavior", {
  tree <- read.tree(text = "((Homo_sapiens_ENSG00001:1,Mus_musculus_ENSMUSG00001:1)n1:1,Danio_rerio_ENSDARG00001:2)root;")
  default_result <- get_species_names(tree)
  legacy_parser <- build_species_parser("legacy")
  explicit_result <- get_species_names(tree, species_parser = legacy_parser)

  expect_equal(default_result, explicit_result)
  expect_equal(default_result, c("Homo_sapiens", "Mus_musculus", "Danio_rerio"))
  expect_no_error(validate_gene_tip_labels(tree$tip.label))
})

test_that("qualified parser accepts qualified species labels", {
  parser <- build_species_parser("qualified")
  species_tips <- c(
    "Dictyostelium_discoideum_cf",
    "Arabidopsis_thaliana_subsp_lyrata",
    "Escherichia_sp_K12"
  )
  gene_tips <- c(
    "Dictyostelium_discoideum_cf_gene123",
    "Arabidopsis_thaliana_subsp_lyrata_gene456",
    "Escherichia_sp_K12_gene789"
  )

  expect_no_error(
    validate_gene_tip_labels(
      gene_tips,
      species_parser = parser,
      species_tree_labels = species_tips
    )
  )
  expect_equal(
    species_parser_get_gene_species(
      parser,
      gene_tips,
      species_tree_labels = species_tips
    ),
    species_tips
  )
  expect_equal(
    leaf2species(
      gene_tips,
      species_parser = parser,
      species_tree_labels = species_tips
    ),
    c(
      "Dictyostelium discoideum cf",
      "Arabidopsis thaliana subsp lyrata",
      "Escherichia sp K12"
    )
  )
})

test_that("qualified_gg remains an alias of qualified", {
  expect_equal(
    build_species_parser("qualified_gg")[["name"]],
    "qualified"
  )
})

test_that("regex parser extracts species labels from custom tip format", {
  parser <- build_species_parser(
    "regex",
    species_regex = "^[^|]+\\|([^|]+)\\|[^|]+$"
  )
  gene_tips <- c(
    "tip1|Dictyostelium_discoideum_cf|gene123",
    "tip2|Homo_sapiens|gene456"
  )

  expect_no_error(validate_gene_tip_labels(gene_tips, species_parser = parser))
  expect_equal(
    species_parser_get_gene_species(parser, gene_tips),
    c("Dictyostelium_discoideum_cf", "Homo_sapiens")
  )
})

test_that("map parser extracts species labels from TSV", {
  map_file <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "label\tspecies\ttaxonomy_query",
      "seq1\tDictyostelium_discoideum_cf\tDictyostelium discoideum cf",
      "seq2\tHomo_sapiens\tHomo sapiens"
    ),
    map_file
  )
  on.exit(unlink(map_file))

  parser <- build_species_parser("map", species_map_tsv = map_file)
  gene_tips <- c("seq1", "seq2")

  expect_no_error(validate_gene_tip_labels(gene_tips, species_parser = parser))
  expect_equal(
    species_parser_get_gene_species(parser, gene_tips),
    c("Dictyostelium_discoideum_cf", "Homo_sapiens")
  )
  expect_equal(
    leaf2species(gene_tips, species_parser = parser),
    c("Dictyostelium discoideum cf", "Homo sapiens")
  )
})

test_that("RADTE allS mode maps mixed legacy and qualified species labels", {
  sp_file <- tempfile(fileext = ".nwk")
  writeLines(
    "((Homo_sapiens:10,Dictyostelium_discoideum_cf:10)n1:10,Arabidopsis_thaliana_subsp_lyrata:20)root;",
    sp_file
  )
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines(
    "((Homo_sapiens_gene1:0.1,Dictyostelium_discoideum_cf_gene123:0.1)n1:0.2,Arabidopsis_thaliana_subsp_lyrata_gene456:0.3)n2;",
    gn_file
  )
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "0\t0\t0\t0\t3\t5\t3\t0.1,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
      "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
      "",
      "#D\tDuplication\tL.Bound\tU.Bound"
    ),
    parsable_file
  )
  on.exit(unlink(parsable_file), add = TRUE)

  result <- run_radte_species_parser(
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--species-parser=qualified",
    "--max_age=1000",
    "--chronos_lambda=1",
    "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  expect_equal(result$exit_code, 0)
  expect_true(inherits(result$output_tree, "phylo"))
  expect_setequal(
    result$output_tree$tip.label,
    c(
      "Homo_sapiens_gene1",
      "Dictyostelium_discoideum_cf_gene123",
      "Arabidopsis_thaliana_subsp_lyrata_gene456"
    )
  )
})
