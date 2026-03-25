test_that("species node bounds TSV parses exact and interval formats", {
  interval_file <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "node\tage_min\tage_max",
      "n1\t5\t7",
      "root\t10\t12"
    ),
    interval_file
  )
  on.exit(unlink(interval_file))

  interval_df <- read_species_node_bounds_tsv(interval_file)
  expect_equal(interval_df$node, c("n1", "root"))
  expect_equal(interval_df$age_min, c(5, 10))
  expect_equal(interval_df$age_max, c(7, 12))

  exact_file <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "node\tage",
      "n1\t6",
      "root\t12"
    ),
    exact_file
  )
  on.exit(unlink(exact_file), add = TRUE)

  exact_df <- read_species_node_bounds_tsv(exact_file)
  expect_equal(exact_df$age_min, c(6, 12))
  expect_equal(exact_df$age_max, c(6, 12))
})

test_that("species node bounds merge rejects intervals that exclude tree ages", {
  sp_node_table <- data.frame(
    node = c("A_sp", "B_sp", "n1", "root"),
    age = c(0, 0, 10, 20),
    age_min = c(0, 0, 10, 20),
    age_max = c(0, 0, 10, 20),
    stringsAsFactors = FALSE
  )
  bounds_df <- data.frame(
    node = "n1",
    age_min = 11,
    age_max = 12,
    stringsAsFactors = FALSE
  )

  expect_error(
    merge_species_node_bounds(sp_node_table, bounds_df),
    "inconsistent with the species-tree branch-length ages"
  )
})

test_that("species constraint metadata marks repeated speciation groups", {
  tree <- read.tree(text = "(((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.1,(A_sp_g2:0.1,B_sp_g2:0.1)n2:0.1)d1:0.1,C_sp_g1:0.2)root;")
  n1_num <- get_node_num_by_name(tree, "n1")
  n2_num <- get_node_num_by_name(tree, "n2")
  d1_num <- get_node_num_by_name(tree, "d1")
  root_num <- get_node_num_by_name(tree, "root")

  gn_node_table <- data.frame(
    event = c("S", "S", "D", "S"),
    gn_node = c("n1", "n2", "d1", "root"),
    gn_node_num = c(n1_num, n2_num, d1_num, root_num),
    lower_sp_node = c("sp_ab", "sp_ab", "sp_ab", "sp_root"),
    upper_sp_node = c("sp_ab", "sp_ab", "sp_root", "sp_root"),
    lower_age = c(5, 5, 5, 10),
    upper_age = c(7, 7, 10, 12),
    stringsAsFactors = FALSE
  )

  annotated <- annotate_species_constraint_groups(gn_node_table, tree)
  expect_equal(annotated$constraint_sp_node[annotated$gn_node == "n1"], "sp_ab")
  expect_equal(annotated$constraint_sp_node[annotated$gn_node == "n2"], "sp_ab")
  expect_equal(sort(unique(na.omit(annotated$shared_speciation_group))), "sp_ab")
  expect_true(all(annotated$shared_speciation_group[annotated$gn_node %in% c("n1", "n2")] == "sp_ab"))
  expect_true(is.na(annotated$shared_speciation_group[annotated$gn_node == "root"]))
})
