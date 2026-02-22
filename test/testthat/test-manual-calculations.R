# Tests with manually-calculated expected values for toy inputs
#
# Each test defines a small input, documents the manual calculation in comments,
# and verifies exact numeric/string values (not just structural properties).

radte_script <- file.path(project_root, "radte.r")

# --- pad_branch_length: exact per-edge values ---

test_that("pad_branch_length: exact per-edge values with one zero-length edge", {
  # Input: ((A:0,B:0.1)n1:0.05,C:0.2)root;  pad_size=0.001
  # Tips: A=1, B=2, C=3. Internal: root=4, n1=5
  # Edge to A has length 0.0 < 0.001 -> set to 0.001
  # All other edges >= 0.001 -> unchanged
  tree <- read.tree(text = "((A:0,B:0.1)n1:0.05,C:0.2)root;")
  result <- pad_branch_length(tree, pad_size = 0.001)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 0.001)
  expect_equal(result$edge.length[b_idx], 0.1)
  expect_equal(result$edge.length[c_idx], 0.2)
  expect_equal(result$edge.length[n1_idx], 0.05)
})

test_that("pad_branch_length: exact values with multiple sub-threshold edges", {
  # Input: ((A:0,B:0.0001)n1:1,C:0.0005)root;  pad_size=0.001
  # A:0 < 0.001 -> 0.001
  # B:0.0001 < 0.001 -> 0.001
  # C:0.0005 < 0.001 -> 0.001
  # n1:1 >= 0.001 -> unchanged
  tree <- read.tree(text = "((A:0,B:0.0001)n1:1,C:0.0005)root;")
  result <- pad_branch_length(tree, pad_size = 0.001)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 0.001)
  expect_equal(result$edge.length[b_idx], 0.001)
  expect_equal(result$edge.length[c_idx], 0.001)
  expect_equal(result$edge.length[n1_idx], 1.0)
})

# --- adjust_branch_length_order: exact scaling factor ---

test_that("adjust_branch_length_order: exact 100x scaling (two rounds)", {
  # Input edges: A:0.00001, B:0.1, n1:0.2, C:0.3
  # Round 1: min=0.00001 < 0.001 -> *10 -> [0.0001, 1, 2, 3]
  # Round 2: min=0.0001 < 0.001 -> *10 -> [0.001, 10, 20, 30]
  # Round 3: min=0.001 >= 0.001 -> break
  tree <- read.tree(text = "((A:0.00001,B:0.1)n1:0.2,C:0.3)root;")
  result <- adjust_branch_length_order(tree, min_bl = 0.001)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 0.001)
  expect_equal(result$edge.length[b_idx], 10.0)
  expect_equal(result$edge.length[c_idx], 30.0)
  expect_equal(result$edge.length[n1_idx], 20.0)
})

test_that("adjust_branch_length_order: exact 10x scaling (single round)", {
  # Input edges: A:0.0005, B:0.5, n1:1.0, C:2.0
  # Round 1: min=0.0005 < 0.001 -> *10 -> [0.005, 5, 10, 20]
  # Round 2: min=0.005 >= 0.001 -> break
  tree <- read.tree(text = "((A:0.0005,B:0.5)n1:1.0,C:2.0)root;")
  result <- adjust_branch_length_order(tree, min_bl = 0.001)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 0.005)
  expect_equal(result$edge.length[b_idx], 5.0)
  expect_equal(result$edge.length[c_idx], 20.0)
  expect_equal(result$edge.length[n1_idx], 10.0)
})

test_that("adjust_branch_length_order: no scaling when min >= threshold", {
  # Input edges: A:0.01, B:0.5, n1:1.0, C:2.0
  # min=0.01 >= 0.001 -> no scaling, all unchanged
  tree <- read.tree(text = "((A:0.01,B:0.5)n1:1.0,C:2.0)root;")
  result <- adjust_branch_length_order(tree, min_bl = 0.001)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 0.01)
  expect_equal(result$edge.length[b_idx], 0.5)
  expect_equal(result$edge.length[c_idx], 2.0)
  expect_equal(result$edge.length[n1_idx], 1.0)
})

# --- pad_short_edges: exact redistribution values ---

test_that("pad_short_edges: exact redistribution from parent edge", {
  # Input: ((A:0.0001,B:1)n1:0.5,C:1.5)root;  threshold=0.001
  # Tips: A=1, B=2, C=3. Internal: root=4, n1=5
  #
  # Short edge: A (0.0001 < 0.001)
  # shift_value = 0.001 - 0.0001 = 0.0009
  # Sister of A is B. Parent edge: root->n1 (0.5 >= 0.001+0.0009=0.0019, OK)
  #
  # After:
  #   edge to A:  0.0001 + 0.0009 = 0.001
  #   edge to B:  1.0 + 0.0009 = 1.0009  (sister gets same shift)
  #   edge to n1: 0.5 - 0.0009 = 0.4991  (parent gives up shift)
  #   edge to C:  1.5 (unchanged)
  #
  # Verify ultrametricity preserved:
  #   root->A: 0.4991 + 0.001 = 0.5001 (same as before: 0.5 + 0.0001)
  #   root->B: 0.4991 + 1.0009 = 1.5 (same as before: 0.5 + 1.0)
  #   root->C: 1.5 (unchanged)
  tree <- read.tree(text = "((A:0.0001,B:1)n1:0.5,C:1.5)root;")
  result <- pad_short_edges(tree, threshold = 0.001)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 0.001)
  expect_equal(result$edge.length[b_idx], 1.0009)
  expect_equal(result$edge.length[n1_idx], 0.4991)
  expect_equal(result$edge.length[c_idx], 1.5)
})

test_that("pad_short_edges: exact redistribution at root (flag_root case)", {
  # Input: ((A:1,B:1)n1:0.0001,C:2)root;  threshold=0.001
  # Short edge: n1 (0.0001 < 0.001)
  # shift_value = 0.001 - 0.0001 = 0.0009
  # Sister of n1 is C. Parent of n1 is root -> flag_root=TRUE
  #
  # After (flag_root: add to both, NO parent subtraction):
  #   edge to n1: 0.0001 + 0.0009 = 0.001
  #   edge to C:  2.0 + 0.0009 = 2.0009
  #   edge to A:  1.0 (unchanged)
  #   edge to B:  1.0 (unchanged)
  tree <- read.tree(text = "((A:1,B:1)n1:0.0001,C:2)root;")
  result <- pad_short_edges(tree, threshold = 0.001)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[n1_idx], 0.001)
  expect_equal(result$edge.length[c_idx], 2.0009)
  expect_equal(result$edge.length[a_idx], 1.0)
  expect_equal(result$edge.length[b_idx], 1.0)
})

test_that("pad_short_edges: external_only=TRUE pads only tip edges", {
  # Input: ((A:0.0001,B:1)n1:0.0001,C:2)root;  threshold=0.001, external_only=TRUE
  # Internal edge n1 (0.0001) is short but NOT a target (external_only=TRUE)
  # External edge A (0.0001) is short and IS a target
  #
  # For A: shift=0.0009, sister=B
  # Walk up: parent of A is n1. parent_edge_length=0.0001 < 0.001+0.0009=0.0019
  # Continue up: parent of n1 is root -> flag_root=TRUE
  #
  # After:
  #   edge to A:  0.0001 + 0.0009 = 0.001
  #   edge to B:  1.0 + 0.0009 = 1.0009  (sister)
  #   edge to n1: 0.0001 (unchanged, not a target)
  #   edge to C:  2.0 (unchanged)
  tree <- read.tree(text = "((A:0.0001,B:1)n1:0.0001,C:2)root;")
  result <- pad_short_edges(tree, threshold = 0.001, external_only = TRUE)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 0.001)
  expect_equal(result$edge.length[b_idx], 1.0009)
  expect_equal(result$edge.length[n1_idx], 0.0001)
  expect_equal(result$edge.length[c_idx], 2.0)
})

# --- force_ultrametric: exact edge lengths preserved ---

test_that("force_ultrametric: exact edge lengths preserved for ultrametric input", {
  # ((A:1,B:1)n1:1,C:2)root; is ultrametric (root->A=2, root->B=2, root->C=2)
  # force_ultrametric should return identical edge lengths
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  result <- force_ultrametric(tree)

  a_idx <- which(result$edge[, 2] == 1)
  b_idx <- which(result$edge[, 2] == 2)
  c_idx <- which(result$edge[, 2] == 3)
  n1_idx <- which(result$edge[, 2] == 5)

  expect_equal(result$edge.length[a_idx], 1.0)
  expect_equal(result$edge.length[b_idx], 1.0)
  expect_equal(result$edge.length[c_idx], 2.0)
  expect_equal(result$edge.length[n1_idx], 1.0)
})

# --- leaf2species: exact string output ---

test_that("leaf2species: exact species from multi-part names", {
  # "Homo_sapiens_BRCA1" -> split ["Homo","sapiens","BRCA1"] -> paste("Homo","sapiens") = "Homo sapiens"
  # "Mus_musculus_Tp53_variant" -> split ["Mus","musculus","Tp53","variant"] -> "Mus musculus"
  # "Danio_rerio_gene1" -> "Danio rerio"
  leaves <- c("Homo_sapiens_BRCA1", "Mus_musculus_Tp53_variant", "Danio_rerio_gene1")
  result <- leaf2species(leaves)
  expect_equal(length(result), 3)
  expect_equal(result[1], "Homo sapiens")
  expect_equal(result[2], "Mus musculus")
  expect_equal(result[3], "Danio rerio")
})

test_that("leaf2species: short names are skipped with warning", {
  # "AB" -> split ["AB"] -> length=1 < 3 -> warning, skipped
  # "X_Y" -> split ["X","Y"] -> length=2 < 3 -> warning, skipped
  # Result: empty character vector
  expect_warning(result <- leaf2species(c("AB", "X_Y")))
  expect_equal(length(result), 0)
})

test_that("leaf2species: mixed valid and invalid names", {
  # "Homo_sapiens_BRCA1" -> "Homo sapiens" (valid, length>=3)
  # "AB" -> skipped with warning (length<3)
  expect_warning(result <- leaf2species(c("Homo_sapiens_BRCA1", "AB")))
  expect_equal(length(result), 1)
  expect_equal(result[1], "Homo sapiens")
})

# --- get_species_names: exact extraction ---

test_that("get_species_names: exact extraction from multi-underscore names", {
  # tip "A_B_C" -> split by '_': ["A","B","C"] -> paste0("A","_","B") = "A_B"
  # tip "X_Y_Z_W" -> split by '_': ["X","Y","Z","W"] -> paste0("X","_","Y") = "X_Y"
  tree <- read.tree(text = "(A_B_C:1,X_Y_Z_W:1)root;")
  result <- get_species_names(tree, sep = "_")
  expect_equal(length(result), 2)
  expect_equal(result[1], "A_B")
  expect_equal(result[2], "X_Y")
})

# --- transfer_node_labels: exact label transfer for 4-tip tree ---

test_that("transfer_node_labels: exact labels for 4-tip tree", {
  # from: ((A,B)FROM_AB,(C,D)FROM_CD)FROM_root
  # to:   ((A,B)TO_AB,(C,D)TO_CD)TO_root
  # Clade {A,B} in to matches {A,B} in from -> TO_AB becomes FROM_AB
  # Clade {C,D} in to matches {C,D} in from -> TO_CD becomes FROM_CD
  # Clade {A,B,C,D} in to matches {A,B,C,D} in from -> TO_root becomes FROM_root
  tree_from <- read.tree(text = "((A:1,B:1)FROM_AB:1,(C:1,D:1)FROM_CD:1)FROM_root;")
  tree_to   <- read.tree(text = "((A:2,B:2)TO_AB:2,(C:2,D:2)TO_CD:2)TO_root;")
  result <- transfer_node_labels(tree_from, tree_to)

  # All three labels should be transferred
  expect_true("FROM_root" %in% result$node.label)
  expect_true("FROM_AB" %in% result$node.label)
  expect_true("FROM_CD" %in% result$node.label)
  # No TO_ labels should remain
  expect_false(any(grepl("^TO_", result$node.label)))
})

# --- read_notung_parsable: cell-by-cell verification ---

test_that("read_notung_parsable: exact cell-by-cell for 3 duplications", {
  # Parsable file with 3 duplication events
  # Line format: #D <gn_node> <lower_sp_node> <upper_sp_node>
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "3\t0\t0\t0\t6\t11\t3\t0.01,2.0\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tdup1\tSpeciesA\tancestor1",
    "#D\tdup2\tSpeciesB\tancestor2",
    "#D\tdup3\tSpeciesC\tancestor3"
  ), tmp)
  on.exit(unlink(tmp))

  result <- read_notung_parsable(tmp, mode = "D")

  # Exact dimensions
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("event", "gn_node", "lower_sp_node", "upper_sp_node"))

  # Row 1: exact cell values
  expect_equal(result$event[1], "D")
  expect_equal(result$gn_node[1], "dup1")
  expect_equal(result$lower_sp_node[1], "SpeciesA")
  expect_equal(result$upper_sp_node[1], "ancestor1")

  # Row 2: exact cell values
  expect_equal(result$event[2], "D")
  expect_equal(result$gn_node[2], "dup2")
  expect_equal(result$lower_sp_node[2], "SpeciesB")
  expect_equal(result$upper_sp_node[2], "ancestor2")

  # Row 3: exact cell values
  expect_equal(result$event[3], "D")
  expect_equal(result$gn_node[3], "dup3")
  expect_equal(result$lower_sp_node[3], "SpeciesC")
  expect_equal(result$upper_sp_node[3], "ancestor3")
})

# --- contains_polytomy: exact table() value check ---

test_that("contains_polytomy: exact max children count", {
  # Binary: each parent has exactly 2 children
  binary <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  expect_equal(max(table(binary$edge[, 1])), 2)
  expect_false(contains_polytomy(binary))

  # Polytomy: root has 3 children
  poly <- read.tree(text = "(A:1,B:1,C:1)root;")
  expect_equal(max(table(poly$edge[, 1])), 3)
  expect_true(contains_polytomy(poly))
})

# --- Pipeline NOTUNG+duplication: manual calculation of all output values ---

test_that("RADTE pipeline: exact sp_node_table, gn_node_table, and calibration values", {
  # Species tree: ((A_sp:10,B_sp:10)n1:20,C_sp:30)root;
  #   node ages: A_sp=0, B_sp=0, C_sp=0, n1=10, root=30
  #   (root depth=30, n1 depth from root=20, age=30-20=10)
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)n1:20,C_sp:30)root;", sp_file)
  on.exit(unlink(sp_file))

  # Gene tree with 1 duplication in A_sp lineage
  # Tips: A_sp_g1, A_sp_g2, B_sp_g1, C_sp_g1
  # n1=duplication, n2=speciation at sp_n1, n3=root speciation at sp_root
  gn_file <- tempfile(fileext = ".nwk")
  writeLines("(((A_sp_g1:0.05,A_sp_g2:0.05)n1:0.05,B_sp_g1:0.1)n2:0.2,C_sp_g1:0.3)n3;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  # Parsable: 1 duplication at n1, between species A_sp (lower) and n1 (upper)
  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "1\t0\t0\t0\t4\t7\t3\t0.05,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound",
    "#D\tn1\tA_sp\tn1"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), paste0("radte_manual_", as.integer(Sys.time())))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000",
    "--chronos_lambda=1",
    "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # --- Verify sp_node_table exact ages ---
  # ages = root_depth - node.depth.edgelength
  # root_depth = 30, so: A_sp=0, B_sp=0, C_sp=0, n1=10, root=30
  sp_tsv <- read.delim(file.path(out_dir, "radte_species_tree.tsv"))
  expect_equal(sp_tsv$age[sp_tsv$node == "A_sp"], 0)
  expect_equal(sp_tsv$age[sp_tsv$node == "B_sp"], 0)
  expect_equal(sp_tsv$age[sp_tsv$node == "C_sp"], 0)
  expect_equal(sp_tsv$age[sp_tsv$node == "root"], 30)
  expect_equal(sp_tsv$age[sp_tsv$node == "n1"], 10)

  # --- Verify gn_node_table exact events and age bounds ---
  gn_tsv <- read.delim(file.path(out_dir, "radte_gene_tree.tsv"))

  # Duplication node n1: event=D, lower_sp=A_sp (age 0), upper_sp=n1 (age 10)
  dup_row <- gn_tsv[gn_tsv$gn_node == "n1", ]
  expect_equal(nrow(dup_row), 1)
  expect_equal(dup_row$event, "D")
  expect_equal(dup_row$lower_sp_node, "A_sp")
  expect_equal(dup_row$upper_sp_node, "n1")
  expect_equal(dup_row$lower_age, 0)
  expect_equal(dup_row$upper_age, 10, tolerance = 1e-6)

  # Speciation node n2: mapped to species n1 (MRCA of A_sp,B_sp in gene subtree)
  # event=S, lower_sp=n1, upper_sp=n1, ages both = 10
  spec_n2 <- gn_tsv[gn_tsv$gn_node == "n2", ]
  expect_equal(nrow(spec_n2), 1)
  expect_equal(spec_n2$event, "S")
  expect_equal(spec_n2$lower_sp_node, "n1")
  expect_equal(spec_n2$upper_sp_node, "n1")
  expect_equal(spec_n2$lower_age, 10, tolerance = 1e-6)
  expect_equal(spec_n2$upper_age, 10, tolerance = 1e-6)

  # Root node n3: mapped to species root
  # event=S(R), lower_sp=root, upper_sp=root, ages both = 30
  root_row <- gn_tsv[gn_tsv$gn_node == "n3", ]
  expect_equal(nrow(root_row), 1)
  expect_equal(root_row$event, "S(R)")
  expect_equal(root_row$lower_sp_node, "root")
  expect_equal(root_row$upper_sp_node, "root")
  expect_equal(root_row$lower_age, 30)
  expect_equal(root_row$upper_age, 30)

  # --- Verify calibration mode and table ---
  # Has duplication -> chronos path -> RS/S/R trials
  # RS should succeed (root + speciation nodes)
  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_equal(cal_nodes[1], "RS")

  cal_used <- read.delim(file.path(out_dir, "radte_calibration_used.tsv"))

  # Root calibration: age.min=30, age.max=30, event=R
  root_cal <- cal_used[cal_used$event == "R", ]
  expect_equal(nrow(root_cal), 1)
  expect_equal(root_cal$age.min, 30)
  expect_equal(root_cal$age.max, 30)

  # Speciation calibration: age.min=10, age.max=10, event=S
  spec_cal <- cal_used[cal_used$event == "S", ]
  expect_equal(nrow(spec_cal), 1)
  expect_equal(spec_cal$age.min, 10, tolerance = 1e-6)
  expect_equal(spec_cal$age.max, 10, tolerance = 1e-6)

  # --- Verify output tree structure ---
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(is.ultrametric(out_tree))
  expect_equal(length(out_tree$tip.label), 4)
  expect_true(all(out_tree$edge.length > 0))
  expect_setequal(out_tree$tip.label,
                  c("A_sp_g1", "A_sp_g2", "B_sp_g1", "C_sp_g1"))
})

# --- Pipeline allS: species tree ages transferred directly ---

test_that("RADTE allS pipeline: exact species ages and tree height", {
  # Species tree: ((A_sp:10,B_sp:10)n1:20,C_sp:30)root;
  # Gene tree: no duplications -> allS mode
  # In allS mode, species tree topology and branch lengths are directly used
  # Expected tree height = species root age = 30
  sp_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp:10,B_sp:10)n1:20,C_sp:30)root;", sp_file)
  on.exit(unlink(sp_file))

  gn_file <- tempfile(fileext = ".nwk")
  writeLines("((A_sp_g1:0.1,B_sp_g1:0.1)n1:0.2,C_sp_g1:0.3)n2;", gn_file)
  on.exit(unlink(gn_file), add = TRUE)

  parsable_file <- tempfile(fileext = ".txt")
  writeLines(c(
    "0\t0\t0\t0\t3\t5\t3\t0.1,0.3\tOff\t1\t1\t1.5\t0.0\t3.0\t1.0",
    "nD\tnCD\tnT\tnL\t|L(G)|\t|G|\t|S|\tminEW,maxEW\tRoots\tCand\tFeas\tcD\tcCD\tcT\tcL",
    "",
    "#D\tDuplication\tL.Bound\tU.Bound"
  ), parsable_file)
  on.exit(unlink(parsable_file), add = TRUE)

  out_dir <- file.path(tempdir(), paste0("radte_allS_manual_", as.integer(Sys.time())))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  cmd <- paste(
    "Rscript", shQuote(radte_script),
    paste0("--species_tree=", shQuote(sp_file)),
    paste0("--gene_tree=", shQuote(gn_file)),
    paste0("--notung_parsable=", shQuote(parsable_file)),
    "--max_age=1000",
    "--chronos_lambda=1",
    "--chronos_model=discrete",
    "--pad_short_edge=0.001"
  )

  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(exit_code, 0)

  # allS mode confirmed
  cal_nodes <- readLines(file.path(out_dir, "radte_calibrated_nodes.txt"))
  expect_equal(cal_nodes[1], "allS")

  # Species tree exact ages
  sp_tsv <- read.delim(file.path(out_dir, "radte_species_tree.tsv"))
  expect_equal(sp_tsv$age[sp_tsv$node == "A_sp"], 0)
  expect_equal(sp_tsv$age[sp_tsv$node == "B_sp"], 0)
  expect_equal(sp_tsv$age[sp_tsv$node == "C_sp"], 0)
  expect_equal(sp_tsv$age[sp_tsv$node == "root"], 30)
  expect_equal(sp_tsv$age[sp_tsv$node == "n1"], 10)

  # Output tree: topology from species tree, gene tip labels
  out_tree <- read.tree(file.path(out_dir, "radte_gene_tree_output.nwk"))
  expect_true(is.ultrametric(out_tree))
  expect_equal(length(out_tree$tip.label), 3)
  expect_setequal(out_tree$tip.label, c("A_sp_g1", "B_sp_g1", "C_sp_g1"))

  # Tree height should equal species tree root age = 30
  tree_height <- max(ape::node.depth.edgelength(out_tree))
  expect_equal(tree_height, 30, tolerance = 0.01)

  # Individual branch lengths should match species tree:
  # root->n1: 20, n1->A_sp_g1: 10, n1->B_sp_g1: 10, root->C_sp_g1: 30
  a_idx <- which(out_tree$edge[, 2] == which(out_tree$tip.label == "A_sp_g1"))
  b_idx <- which(out_tree$edge[, 2] == which(out_tree$tip.label == "B_sp_g1"))
  c_idx <- which(out_tree$edge[, 2] == which(out_tree$tip.label == "C_sp_g1"))
  expect_equal(out_tree$edge.length[a_idx], 10, tolerance = 0.01)
  expect_equal(out_tree$edge.length[b_idx], 10, tolerance = 0.01)
  expect_equal(out_tree$edge.length[c_idx], 30, tolerance = 0.01)
})
