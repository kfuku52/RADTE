# Test tree: ((A:1,B:1)n1:1,C:2)root;
# Tips: A=1, B=2, C=3; Internal: root=4, n1=5

make_test_tree <- function() {
  read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
}

test_that("get_root_num returns the root node number", {
  tree <- make_test_tree()
  root <- get_root_num(tree)
  expect_equal(root, 4)  # 3 tips + 1 = 4
})

test_that("get_node_name_by_num maps tip nodes correctly", {
  tree <- make_test_tree()
  expect_equal(get_node_name_by_num(tree, 1), "A")
  expect_equal(get_node_name_by_num(tree, 2), "B")
  expect_equal(get_node_name_by_num(tree, 3), "C")
})

test_that("get_node_name_by_num maps internal nodes correctly", {
  tree <- make_test_tree()
  expect_equal(get_node_name_by_num(tree, 4), "root")
  expect_equal(get_node_name_by_num(tree, 5), "n1")
})

test_that("get_node_num_by_name maps tip names correctly", {
  tree <- make_test_tree()
  expect_equal(get_node_num_by_name(tree, "A"), 1)
  expect_equal(get_node_num_by_name(tree, "B"), 2)
  expect_equal(get_node_num_by_name(tree, "C"), 3)
})

test_that("get_node_num_by_name maps internal node names correctly", {
  tree <- make_test_tree()
  expect_equal(get_node_num_by_name(tree, "root"), 4)
  expect_equal(get_node_num_by_name(tree, "n1"), 5)
})

test_that("get_parent_num returns correct parent for tips", {
  tree <- make_test_tree()
  # A and B are children of n1 (node 5)
  expect_equal(get_parent_num(tree, 1), 5)  # parent of A is n1
  expect_equal(get_parent_num(tree, 2), 5)  # parent of B is n1
  # C is child of root (node 4)
  expect_equal(get_parent_num(tree, 3), 4)  # parent of C is root
})

test_that("get_parent_num returns correct parent for internal nodes", {
  tree <- make_test_tree()
  # n1 is child of root
  expect_equal(get_parent_num(tree, 5), 4)
})

test_that("get_sister_num returns correct sister for tips", {
  tree <- make_test_tree()
  # A's sister is B
  expect_equal(get_sister_num(tree, 1), 2)
  # B's sister is A
  expect_equal(get_sister_num(tree, 2), 1)
  # C's sister is n1
  expect_equal(get_sister_num(tree, 3), 5)
})

test_that("get_sister_num returns correct sister for internal nodes", {
  tree <- make_test_tree()
  # n1's sister is C
  expect_equal(get_sister_num(tree, 5), 3)
})

test_that("get_ancestor_num returns all ancestors to root", {
  tree <- make_test_tree()
  # Ancestors of A (node 1): n1 (5), root (4)
  ancestors <- get_ancestor_num(tree, 1)
  expect_equal(ancestors, c(5, 4))
})

test_that("get_ancestor_num returns empty for root", {
  tree <- make_test_tree()
  ancestors <- get_ancestor_num(tree, 4)
  expect_equal(length(ancestors), 0)
})

test_that("get_ancestor_num returns parent for direct child of root", {
  tree <- make_test_tree()
  # Ancestors of C (node 3): root (4)
  ancestors <- get_ancestor_num(tree, 3)
  expect_equal(ancestors, c(4))
})
