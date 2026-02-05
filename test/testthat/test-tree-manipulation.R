# --- contains_polytomy ---

test_that("contains_polytomy returns FALSE for binary tree", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  expect_false(contains_polytomy(tree))
})

test_that("contains_polytomy returns TRUE for multifurcating tree", {
  tree <- read.tree(text = "(A:1,B:1,C:1)root;")
  expect_true(contains_polytomy(tree))
})

# --- pad_branch_length ---

test_that("pad_branch_length pads zero-length branches", {
  tree <- read.tree(text = "((A:0,B:1)n1:1,C:2)root;")
  result <- pad_branch_length(tree, pad_size = 1e-6)
  expect_true(all(result$edge.length >= 1e-6))
})

test_that("pad_branch_length does not modify branches above threshold", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  result <- pad_branch_length(tree, pad_size = 1e-6)
  expect_equal(tree$edge.length, result$edge.length)
})

test_that("pad_branch_length pads with custom pad_size", {
  tree <- read.tree(text = "((A:0.0001,B:1)n1:1,C:2)root;")
  result <- pad_branch_length(tree, pad_size = 0.001)
  min_bl <- min(result$edge.length)
  expect_gte(min_bl, 0.001)
})

# --- adjust_branch_length_order ---

test_that("adjust_branch_length_order scales up small branch lengths", {
  tree <- read.tree(text = "((A:1e-8,B:1e-7)n1:1e-7,C:1e-6)root;")
  result <- adjust_branch_length_order(tree, min_bl = 1e-6)
  expect_gte(min(result$edge.length), 1e-6)
})

test_that("adjust_branch_length_order does not modify adequate branch lengths", {
  tree <- read.tree(text = "((A:1,B:2)n1:1,C:3)root;")
  result <- adjust_branch_length_order(tree, min_bl = 1e-6)
  expect_equal(tree$edge.length, result$edge.length)
})

test_that("adjust_branch_length_order errors on zero branch length", {
  tree <- read.tree(text = "((A:0,B:1)n1:1,C:2)root;")
  expect_error(adjust_branch_length_order(tree, min_bl = 1e-6))
})

# --- force_ultrametric ---

test_that("force_ultrametric returns ultrametric tree unchanged", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  result <- force_ultrametric(tree)
  expect_true(is.ultrametric(result))
  expect_equal(tree$edge.length, result$edge.length)
})

test_that("force_ultrametric adjusts non-ultrametric tree", {
  tree <- read.tree(text = "((A:1,B:1.001)n1:1,C:2)root;")
  result <- force_ultrametric(tree)
  expect_true(is.ultrametric(result))
})

# --- pad_short_edges ---

test_that("pad_short_edges extends short edges", {
  tree <- read.tree(text = "((A:0.0001,B:1)n1:1,C:2)root;")
  result <- pad_short_edges(tree, threshold = 0.001)
  expect_gte(min(result$edge.length), 0.001)
})

test_that("pad_short_edges does not modify tree without short edges", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  result <- pad_short_edges(tree, threshold = 1e-6)
  expect_equal(tree$edge.length, result$edge.length)
})

test_that("pad_short_edges with external_only=TRUE only pads external edges", {
  tree <- read.tree(text = "((A:0.0001,B:1)n1:1,C:2)root;")
  result <- pad_short_edges(tree, threshold = 0.001, external_only = TRUE)
  # The short external edge (A) should be padded
  a_idx <- which(tree$edge[, 2] == 1)
  expect_gte(result$edge.length[a_idx], 0.001)
})
