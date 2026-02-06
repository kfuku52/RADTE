# Edge case tests derived from GitHub issues #4, #5, #6, #8

# --- read_notung_parsable: no duplication events ---

test_that("read_notung_parsable returns empty df when no duplication events exist", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t0\t0\t0\t0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), tmp)
  on.exit(unlink(tmp))

  result <- read_notung_parsable(tmp, mode = "D")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

# --- read_notung_parsable: tab-delimited fields (issue #3 data format) ---

test_that("read_notung_parsable handles tab-delimited parsable files", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t3\t5\t3\t0.1,1.0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\tSpeciesA\troot"
  ), tmp)
  on.exit(unlink(tmp))

  result <- read_notung_parsable(tmp, mode = "D")
  expect_equal(nrow(result), 1)
  expect_equal(result$gn_node[1], "n1")
  expect_equal(result$lower_sp_node[1], "SpeciesA")
  expect_equal(result$upper_sp_node[1], "root")
})

# --- force_ultrametric: large adjustment exceeds threshold ---

test_that("force_ultrametric errors when adjustment exceeds threshold (issue #6 related)", {
  # Create a highly non-ultrametric tree
  tree <- read.tree(text = "((A:1,B:5)n1:1,C:2)root;")
  expect_error(force_ultrametric(tree, stop_if_larger_change = 0.01))
})

# --- pad_short_edges: short edge at root's direct child ---

test_that("pad_short_edges handles short edge at root child (issue #5 related)", {
  tree <- read.tree(text = "((A:1,B:1)n1:0.0001,C:2)root;")
  result <- pad_short_edges(tree, threshold = 0.001)
  # n1 edge should now be at least threshold
  n1_idx <- which(tree$edge[, 2] == get_node_num_by_name(tree, "n1"))
  expect_gte(result$edge.length[n1_idx], 0.001)
})

# --- contains_polytomy: larger trees ---

test_that("contains_polytomy detects polytomy in larger tree", {
  tree <- read.tree(text = "((A:1,B:1,C:1)n1:1,(D:1,E:1)n2:1)root;")
  expect_true(contains_polytomy(tree))
})

test_that("contains_polytomy returns FALSE for larger binary tree", {
  tree <- read.tree(text = "(((A:1,B:1)n1:1,C:1)n2:1,(D:1,E:1)n3:1)root;")
  expect_false(contains_polytomy(tree))
})

# --- pad_branch_length: multiple zero-length branches ---

test_that("pad_branch_length pads multiple zero-length branches (issue #5 related)", {
  tree <- read.tree(text = "((A:0,B:0)n1:0.1,C:0)root;")
  result <- pad_branch_length(tree, pad_size = 0.001)
  expect_true(all(result$edge.length >= 0.001))
  # The non-zero branch should remain unchanged
  n1_idx <- which(result$edge[, 2] == get_node_num_by_name(result, "n1"))
  expect_equal(result$edge.length[n1_idx], 0.1)
})

# --- adjust_branch_length_order: negative branch length ---

test_that("adjust_branch_length_order errors on negative branch lengths", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  tree$edge.length[1] <- -0.1
  expect_error(adjust_branch_length_order(tree, min_bl = 1e-6))
})

# --- leaf2species: various name formats ---

test_that("leaf2species handles names with more than 3 underscore-separated parts", {
  leaves <- c("Homo_sapiens_ENSG00000_102144")
  result <- leaf2species(leaves)
  expect_equal(result, "Homo sapiens")
})

# --- get_species_names: single-tip tree ---

test_that("get_species_names works with two-tip tree", {
  tree <- read.tree(text = "(Homo_sapiens_ENSG00001:1,Mus_musculus_ENSMUSG00001:1)root;")
  result <- get_species_names(tree)
  expect_equal(result, c("Homo_sapiens", "Mus_musculus"))
})

# --- transfer_node_labels: different leaf sets ---

test_that("transfer_node_labels keeps original labels for non-matching clades", {
  tree_from <- read.tree(text = "((X:1,Y:1)FROM_n1:1,Z:2)FROM_root;")
  tree_to   <- read.tree(text = "((A:1,B:1)TO_n1:1,C:2)TO_root;")
  result <- transfer_node_labels(tree_from, tree_to)
  # No clades match, so labels should remain unchanged
  expect_equal(result$node.label, c("TO_root", "TO_n1"))
})

# --- save_tree_pdf: output file creation ---

test_that("save_tree_pdf creates a valid PDF file", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  pdf_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(pdf_file))

  save_tree_pdf(tree, file = pdf_file)
  expect_true(file.exists(pdf_file))
  expect_gt(file.info(pdf_file)$size, 0)
})

test_that("save_tree_pdf creates PDF with edge colors and ages", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  pdf_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(pdf_file))

  ec <- list("blue" = c(4, 5), "red" = c(1, 2, 3))
  save_tree_pdf(tree, file = pdf_file, show.age = TRUE, edge_colors = ec)
  expect_true(file.exists(pdf_file))
  expect_gt(file.info(pdf_file)$size, 0)
})
