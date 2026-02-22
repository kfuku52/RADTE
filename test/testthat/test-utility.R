# --- get_species_names ---

test_that("get_species_names extracts genus_species from tip labels", {
  tree <- read.tree(text = "((Homo_sapiens_ENSG00001:1,Mus_musculus_ENSMUSG00001:1)n1:1,Danio_rerio_ENSDARG00001:2)root;")
  result <- get_species_names(tree)
  expect_equal(result, c("Homo_sapiens", "Mus_musculus", "Danio_rerio"))
})

test_that("get_species_names handles multiple underscores in gene ID", {
  tree <- read.tree(text = "((Homo_sapiens_ENSG_00001:1,Mus_musculus_ENSMUSG_00001:1)n1:1,Danio_rerio_ENSDARG_00001:2)root;")
  result <- get_species_names(tree)
  expect_equal(result, c("Homo_sapiens", "Mus_musculus", "Danio_rerio"))
})

# --- leaf2species ---

test_that("leaf2species converts leaf names to species names", {
  leaves <- c("Homo_sapiens_ENSG00001", "Mus_musculus_ENSMUSG00001")
  result <- leaf2species(leaves)
  expect_equal(result, c("Homo sapiens", "Mus musculus"))
})

test_that("leaf2species warns for names with fewer than 3 parts", {
  leaves <- c("SingleName")
  expect_warning(leaf2species(leaves))
})

test_that("leaf2species handles empty input", {
  result <- leaf2species(character(0))
  expect_equal(length(result), 0)
})

# --- transfer_node_labels ---

test_that("transfer_node_labels transfers labels between trees with same topology", {
  tree_from <- read.tree(text = "((A:1,B:1)FROM_n1:1,C:2)FROM_root;")
  tree_to   <- read.tree(text = "((A:2,B:2)TO_n1:2,C:3)TO_root;")
  result <- transfer_node_labels(tree_from, tree_to)
  expect_equal(result$node.label, c("FROM_root", "FROM_n1"))
})

test_that("transfer_node_labels preserves labels when no match found", {
  tree_from <- read.tree(text = "((A:1,B:1)FROM_n1:1,D:2)FROM_root;")
  tree_to   <- read.tree(text = "((A:2,B:2)TO_n1:2,C:3)TO_root;")
  result <- transfer_node_labels(tree_from, tree_to)
  # n1 clade {A,B} matches, so it transfers; root clade differs, stays
  expect_equal(result$node.label[2], "FROM_n1")
})

# --- check_gn_node_name_uniqueness ---

test_that("check_gn_node_name_uniqueness passes for unique names", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  gn_table <- data.frame(gn_node = c("n1"), stringsAsFactors = FALSE)
  # Should not error
  expect_no_error(check_gn_node_name_uniqueness(gn_table, tree))
})

test_that("check_gn_node_name_uniqueness errors for non-existent node names", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  gn_table <- data.frame(gn_node = c("nonexistent"), stringsAsFactors = FALSE)
  expect_error(check_gn_node_name_uniqueness(gn_table, tree))
})

# --- ensure_root_event_tag ---

test_that("ensure_root_event_tag appends root marker only when missing", {
  expect_equal(ensure_root_event_tag("S"), "S(R)")
  expect_equal(ensure_root_event_tag("D"), "D(R)")
  expect_equal(ensure_root_event_tag("S(R)"), "S(R)")
  expect_equal(ensure_root_event_tag("D(R)"), "D(R)")
  expect_equal(ensure_root_event_tag("R"), "R")
})
