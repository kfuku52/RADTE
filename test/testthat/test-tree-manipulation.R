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

# --- normalize_edge_length_range ---

test_that("normalize_edge_length_range scales down large edge lengths", {
  tree <- read.tree(text = "((A:1e-15,B:1e10)n1:1,(C:1,D:1)n2:1)root;")
  result <- normalize_edge_length_range(tree, max_edge = 1000, min_edge = 1e-8)
  expect_lte(max(result$edge.length), 1000)
  expect_gte(min(result$edge.length), 1e-8)
})

test_that("normalize_edge_length_range does not modify normal trees", {
  tree <- read.tree(text = "((A:1,B:2)n1:1,C:3)root;")
  result <- normalize_edge_length_range(tree, max_edge = 1000, min_edge = 1e-8)
  expect_equal(tree$edge.length, result$edge.length)
})

test_that("normalize_edge_length_range preserves proportions", {
  tree <- read.tree(text = "((A:1000,B:2000)n1:3000,C:4000)root;")
  result <- normalize_edge_length_range(tree, max_edge = 100, min_edge = 1e-8)
  # Original proportions: 1:2:3:4, should be preserved
  expect_equal(max(result$edge.length), 100)
  expect_equal(result$edge.length[result$edge.length == 100] /
               result$edge.length[result$edge.length == min(result$edge.length)],
               4000 / 1000)
})

test_that("normalize_edge_length_range enables chronos on extreme trees", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1e-15,B:1e10)n1:1,(C:1,D:1)n2:1)root;")
  tree_norm <- normalize_edge_length_range(tree, max_edge = 1000, min_edge = 1e-8)
  calibration <- makeChronosCalib(tree_norm, node = "root", age.min = 100, age.max = 100)
  # Should not error
  expect_no_error(suppressMessages(suppressWarnings(
    chronos(tree_norm, lambda = 1, model = "discrete", calibration = calibration)
  )))
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
