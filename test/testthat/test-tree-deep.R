# Tests with deeper/larger trees for navigation functions

# Test tree: (((A:1,B:1)n3:1,(C:1,D:1)n4:1)n2:1,E:3)n1;
# Tips: A=1,B=2,C=3,D=4,E=5; Internal: n1=6,n2=7,n3=8,n4=9

make_deep_tree <- function() {
  read.tree(text = "(((A:1,B:1)n3:1,(C:1,D:1)n4:1)n2:1,E:3)n1;")
}

test_that("get_root_num on 5-tip tree", {
  tree <- make_deep_tree()
  expect_equal(get_root_num(tree), 6)
})

test_that("get_ancestor_num returns full chain from deepest tip", {
  tree <- make_deep_tree()
  # A (node 1) -> n3 (8) -> n2 (7) -> n1 (6)
  ancestors <- get_ancestor_num(tree, 1)
  expect_equal(ancestors, c(8, 7, 6))
})

test_that("get_ancestor_num returns shorter chain from shallow tip", {
  tree <- make_deep_tree()
  # E (node 5) -> n1 (6)
  ancestors <- get_ancestor_num(tree, 5)
  expect_equal(ancestors, c(6))
})

test_that("get_sister_num for deep internal node", {
  tree <- make_deep_tree()
  # n3 (8) and n4 (9) are sisters
  expect_equal(get_sister_num(tree, 8), 9)
  expect_equal(get_sister_num(tree, 9), 8)
})

test_that("get_parent_num for deep internal node", {
  tree <- make_deep_tree()
  # parent of n3 (8) is n2 (7)
  expect_equal(get_parent_num(tree, 8), 7)
  # parent of n2 (7) is n1 (6)
  expect_equal(get_parent_num(tree, 7), 6)
})

test_that("get_node_num_by_name returns empty for non-existent name", {
  tree <- make_deep_tree()
  result <- get_node_num_by_name(tree, "nonexistent")
  expect_equal(length(result), 0)
})

test_that("get_node_name_by_num and get_node_num_by_name are inverse", {
  tree <- make_deep_tree()
  for (i in 1:9) {
    name <- get_node_name_by_num(tree, i)
    num <- get_node_num_by_name(tree, name)
    expect_equal(num, i)
  }
})

# --- pad_short_edges with cascading short edges ---

test_that("pad_short_edges handles cascading short edges", {
  # Both n3->A and n2->n3 are short
  tree <- read.tree(text = "(((A:0.0001,B:1)n3:0.0001,(C:1,D:1)n4:1)n2:1,E:3)n1;")
  result <- pad_short_edges(tree, threshold = 0.001)
  # All edges should be >= threshold now
  a_idx <- which(result$edge[, 2] == 1)
  expect_gte(result$edge.length[a_idx], 0.001)
})

test_that("pad_short_edges flag_root case: short edge at root child", {
  # n2 is root's child and has short edge, but root has no parent
  tree <- read.tree(text = "((A:1,B:1)n2:0.0001,C:2)n1;")
  result <- pad_short_edges(tree, threshold = 0.001)
  n2_idx <- which(result$edge[, 2] == get_node_num_by_name(result, "n2"))
  expect_gte(result$edge.length[n2_idx], 0.001)
})

# --- contains_polytomy on deeper trees ---

test_that("contains_polytomy detects polytomy at non-root level", {
  tree <- read.tree(text = "((A:1,B:1,C:1)n2:1,D:2)n1;")
  expect_true(contains_polytomy(tree))
})

test_that("contains_polytomy on balanced binary tree", {
  tree <- make_deep_tree()
  expect_false(contains_polytomy(tree))
})

# --- force_ultrametric on deeper tree ---

test_that("force_ultrametric works on deeper ultrametric tree", {
  tree <- make_deep_tree()
  result <- force_ultrametric(tree)
  expect_true(is.ultrametric(result))
})
