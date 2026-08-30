test_that("deep-tree matching and constraint propagation preserve feasible bounds", {
  skip_if(
    tolower(Sys.getenv("RADTE_SKIP_PERFORMANCE_TESTS", "false")) %in% c("1", "true", "yes"),
    "RADTE_SKIP_PERFORMANCE_TESTS requested"
  )
  tree <- stree(1500, type = "left")
  internal <- seq.int(Ntip(tree) + 1L, Ntip(tree) + tree$Nnode)
  constraints <- data.frame(
    gn_node_num = internal,
    lower_age = rev(seq_along(internal)),
    upper_age = rev(seq_along(internal)) + 0.5,
    stringsAsFactors = FALSE
  )
  root <- get_root_num(tree)

  {
    matches <- match_internal_nodes_by_clade(tree, tree)
    stabilized <- stabilize_descendant_constraints(constraints, tree, root)
  }

  expect_false(anyNA(matches))
  expect_equal(nrow(stabilized$gn_node_table), nrow(constraints))
  expect_identical(stabilized$gn_node_table, constraints)
})
