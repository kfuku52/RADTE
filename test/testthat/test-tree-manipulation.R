# --- contains_polytomy ---

test_that("contains_polytomy returns FALSE for binary tree", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  expect_false(contains_polytomy(tree))
})

test_that("contains_polytomy returns TRUE for multifurcating tree", {
  tree <- read.tree(text = "(A:1,B:1,C:1)root;")
  expect_true(contains_polytomy(tree))
})

test_that("contains_polytomy returns TRUE for tree with singleton internal node", {
  tree <- read.tree(text = "((A:1)n1:1,B:1)root;")
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

# --- chronos constraint stabilization ---

test_that("find_descendant_constraint_conflicts detects same-age descendant constraints", {
  tree <- read.tree(text = "(((A:0.1,B:0.1)n1:0.1,C:0.2)n2:0.3,D:0.5)root;")
  root_num <- get_root_num(tree)
  node_map <- setNames((length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), tree$node.label)

  gn_node_table <- data.frame(
    event = c("S(R)", "S", "S"),
    gn_node = c("root", "n2", "n1"),
    gn_node_num = c(root_num, node_map[["n2"]], node_map[["n1"]]),
    lower_age = c(100, 80, 80),
    upper_age = c(100, 80, 80),
    stringsAsFactors = FALSE
  )

  conflicts <- find_descendant_constraint_conflicts(gn_node_table, tree, root_num)
  expect_equal(nrow(conflicts), 1)
  expect_equal(conflicts$node[1], node_map[["n1"]])
})

test_that("stabilize_descendant_constraints rejects impossible bounds without changing them", {
  tree <- read.tree(text = "(((A:0.1,B:0.1)n1:0.1,C:0.2)n2:0.3,D:0.5)root;")
  root_num <- get_root_num(tree)
  node_map <- setNames((length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), tree$node.label)

  gn_node_table <- data.frame(
    event = c("S(R)", "S", "S"),
    gn_node = c("root", "n2", "n1"),
    gn_node_num = c(root_num, node_map[["n2"]], node_map[["n1"]]),
    lower_age = c(100, 80, 80),
    upper_age = c(100, 80, 80),
    stringsAsFactors = FALSE
  )

  original <- gn_node_table
  expect_error(stabilize_descendant_constraints(gn_node_table, tree, root_num), "temporally infeasible")
  expect_identical(gn_node_table, original)
})

test_that("detect_chronos_failure_risks flags tight constraints and extreme edge ranges", {
  tree <- read.tree(text = "((A:1e-10,B:1e5)n1:0.1,C:0.2)root;")
  root_num <- get_root_num(tree)
  node_map <- setNames((length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), tree$node.label)

  calibration <- data.frame(
    node = as.integer(c(root_num, node_map[["n1"]])),
    age.min = c(100, 50),
    age.max = c(100, 50),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  risks <- detect_chronos_failure_risks(tree, calibration, root_num)
  expect_true(isTRUE(risks$risk_flags[["extreme_edge_ratio"]]))
  expect_true(isTRUE(risks$risk_flags[["tight_nonroot_constraints"]]))
  expect_equal(length(risks$tight_nodes), 1)
})

test_that("make_chronos_control_profiles uses a fast profile before high fallback", {
  profiles <- make_chronos_control_profiles()

  expect_equal(names(profiles), c("fast", "high-fallback"))
  expect_equal(profiles$fast$iter.max, 250)
  expect_equal(profiles$fast$eval.max, 250)
  expect_equal(profiles$fast$dual.iter.max, 20)
  expect_equal(profiles[["high-fallback"]]$iter.max, 100000)
  expect_equal(profiles[["high-fallback"]]$eval.max, 100000)
  expect_equal(profiles[["high-fallback"]]$dual.iter.max, 200)
})

test_that("make_chronos_control_profiles can disable high fallback", {
  profiles <- make_chronos_control_profiles(
    iter_max = 500,
    eval_max = 600,
    dual_iter_max = 10,
    enable_high_fallback = FALSE
  )

  expect_equal(names(profiles), "fast")
  expect_equal(profiles$fast$iter.max, 500)
  expect_equal(profiles$fast$eval.max, 600)
  expect_equal(profiles$fast$dual.iter.max, 10)
})


test_that("validate_duplication_nodes_internal returns invisibly when no duplication rows exist", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  gn_node_table <- data.frame(
    event = c("S(R)", "S"),
    gn_node = c("root", "n1"),
    gn_node_num = c(get_root_num(tree), get_node_num_by_name(tree, "n1")),
    stringsAsFactors = FALSE
  )
  call_result <- withVisible(validate_duplication_nodes_internal(gn_node_table, tree))
  expect_false(call_result$visible)
  expect_null(call_result$value)
})

test_that("validate_gn_node_table returns invisibly for empty table", {
  empty_table <- data.frame(
    gn_node = character(),
    gn_node_num = integer(),
    stringsAsFactors = FALSE
  )
  call_result <- withVisible(validate_gn_node_table(empty_table))
  expect_false(call_result$visible)
  expect_null(call_result$value)
})

test_that("run_chronos_with_restarts succeeds on a valid calibration", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  root_num <- get_root_num(tree)
  calibration <- data.frame(
    node = as.integer(root_num),
    age.min = 50,
    age.max = 50,
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  ctl <- chronos.control()
  ctl$iter.max <- 5000
  ctl$eval.max <- 5000
  ctl$dual.iter.max <- 50
  out <- suppressWarnings(
    run_chronos_with_restarts(
      phy = tree,
      calibration = calibration,
      chronos_control = ctl,
      chronos_lambda = 1,
      chronos_model = "discrete",
      context_label = "test-valid",
      max_restarts = 2,
      seed_base = 1
    )
  )
  expect_true(out$success)
  expect_true(inherits(out$chronos_out, "chronos"))
  expect_equal(out$used_seed, 1)
})

test_that("run_chronos_with_restarts reports failure on contradictory calibrations", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  calibration <- data.frame(
    node = as.integer(c(root_num, n1_num)),
    age.min = c(10, 20),
    age.max = c(10, 20),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  ctl <- chronos.control()
  ctl$iter.max <- 5000
  ctl$eval.max <- 5000
  ctl$dual.iter.max <- 50
  out <- run_chronos_with_restarts(
    phy = tree,
    calibration = calibration,
    chronos_control = ctl,
    chronos_lambda = 1,
    chronos_model = "discrete",
    context_label = "test-invalid",
    max_restarts = 1,
    seed_base = 1
  )
  expect_false(out$success)
  expect_true("try-error" %in% class(out$chronos_out))
})

test_that("run_chronos_with_restarts converts elapsed-time-limit errors into timeout failures", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  calibration <- data.frame(
    node = as.integer(get_root_num(tree)),
    age.min = 10,
    age.max = 10,
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  had_chronos <- exists("radte_chronos", envir = environment(run_chronos_with_restarts), inherits = FALSE)
  old_chronos <- NULL
  if (had_chronos) {
    old_chronos <- get("radte_chronos", envir = environment(run_chronos_with_restarts), inherits = FALSE)
  }
  assign(
    "radte_chronos",
    function(...) {
      stop("reached elapsed time limit")
    },
    envir = environment(run_chronos_with_restarts)
  )
  on.exit({
    if (had_chronos) {
      assign("radte_chronos", old_chronos, envir = environment(run_chronos_with_restarts))
    } else if (exists("radte_chronos", envir = environment(run_chronos_with_restarts), inherits = FALSE)) {
      rm("radte_chronos", envir = environment(run_chronos_with_restarts))
    }
  }, add = TRUE)

  out <- run_chronos_with_restarts(
    phy = tree,
    calibration = calibration,
    chronos_control = list(),
    chronos_lambda = 1,
    chronos_model = "discrete",
    context_label = "test-timeout",
    max_restarts = 1,
    seed_base = 1,
    attempt_timeout_sec = 0.01
  )

  expect_false(out$success)
  expect_true("try-error" %in% class(out$chronos_out))
  expect_match(as.character(out$chronos_out), "timed out", ignore.case = TRUE)
})

test_that("run_chronos_with_restarts respects total chronos time budget", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  calibration <- data.frame(
    node = as.integer(get_root_num(tree)),
    age.min = 10,
    age.max = 10,
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  had_chronos <- exists("radte_chronos", envir = environment(run_chronos_with_restarts), inherits = FALSE)
  old_chronos <- NULL
  if (had_chronos) {
    old_chronos <- get("radte_chronos", envir = environment(run_chronos_with_restarts), inherits = FALSE)
  }
  assign(
    "radte_chronos",
    function(...) {
      stop("chronos should not run when time budget is exhausted")
    },
    envir = environment(run_chronos_with_restarts)
  )
  on.exit({
    if (had_chronos) {
      assign("radte_chronos", old_chronos, envir = environment(run_chronos_with_restarts))
    } else if (exists("radte_chronos", envir = environment(run_chronos_with_restarts), inherits = FALSE)) {
      rm("radte_chronos", envir = environment(run_chronos_with_restarts))
    }
  }, add = TRUE)

  budget <- create_chronos_time_budget(0.01)
  Sys.sleep(0.03)
  out <- run_chronos_with_restarts(
    phy = tree,
    calibration = calibration,
    chronos_control = list(),
    chronos_lambda = 1,
    chronos_model = "discrete",
    context_label = "test-budget",
    max_restarts = 2,
    seed_base = 1,
    attempt_timeout_sec = 1,
    time_budget = budget
  )

  expect_false(out$success)
  expect_true("try-error" %in% class(out$chronos_out))
  expect_match(as.character(out$chronos_out), "time budget", ignore.case = TRUE)
})

test_that("run_chronos_with_restarts bounds an active attempt by the total budget", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  calibration <- data.frame(
    node = as.integer(get_root_num(tree)),
    age.min = 10,
    age.max = 10,
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )

  had_runner <- exists("run_chronos_once", envir = environment(run_chronos_with_restarts), inherits = FALSE)
  old_runner <- NULL
  if (had_runner) {
    old_runner <- get("run_chronos_once", envir = environment(run_chronos_with_restarts), inherits = FALSE)
  }
  observed_timeout <- NA_real_
  assign(
    "run_chronos_once",
    function(..., timeout_sec = Inf) {
      observed_timeout <<- timeout_sec
      structure("mock failure", class = "try-error")
    },
    envir = environment(run_chronos_with_restarts)
  )
  on.exit({
    if (had_runner) {
      assign("run_chronos_once", old_runner, envir = environment(run_chronos_with_restarts))
    } else if (exists("run_chronos_once", envir = environment(run_chronos_with_restarts), inherits = FALSE)) {
      rm("run_chronos_once", envir = environment(run_chronos_with_restarts))
    }
  }, add = TRUE)

  budget <- create_chronos_time_budget(0.5)
  out <- run_chronos_with_restarts(
    phy = tree,
    calibration = calibration,
    chronos_control = list(),
    chronos_lambda = 1,
    chronos_model = "discrete",
    context_label = "test-active-budget",
    max_restarts = 1,
    seed_base = 1,
    attempt_timeout_sec = Inf,
    time_budget = budget
  )

  expect_false(out$success)
  expect_true(is.finite(observed_timeout))
  expect_gt(observed_timeout, 0)
  expect_lte(observed_timeout, 0.5)
})

test_that("run_chronos_retry_pipeline skips empty calibration tables", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  out <- run_chronos_retry_pipeline(
    phy = tree,
    calibration_table = data.frame(
      node = integer(),
      age.min = numeric(),
      age.max = numeric(),
      soft.bounds = logical(),
      stringsAsFactors = FALSE
    ),
    root_num = get_root_num(tree),
    chronos_control = list(),
    chronos_lambda = 1,
    chronos_model = "discrete",
    alternate_strategies = list(list(model = "discrete", lambda = 1, label = "requested")),
    calibration_label = "S"
  )

  expect_false(out$success)
  expect_true(isTRUE(out$skipped))
  expect_true("try-error" %in% class(out$chronos_out))
  expect_match(as.character(out$chronos_out), "No calibration nodes", ignore.case = TRUE)
})

test_that("retry strategies preserve bounds when switching model", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  calibration <- data.frame(node = c(4L, 5L), age.min = c(10, 5), age.max = c(10, 5), soft.bounds = FALSE)
  env <- new.env(parent = environment(run_chronos_retry_pipeline))
  seen <- list()
  env$run_chronos_with_restarts <- function(...) {
    args <- list(...)
    seen[[length(seen) + 1L]] <<- args$calibration
    success <- identical(args$chronos_model, "relaxed")
    list(chronos_out = if (success) tree else structure("retry", class = "try-error"),
         success = success, used_model = args$chronos_model, used_lambda = args$chronos_lambda,
         used_seed = args$seed_base)
  }
  runner <- run_chronos_retry_pipeline
  environment(runner) <- env
  out <- runner(tree, calibration, 4L, list(), 1, "discrete",
                list(list(model = "relaxed", lambda = 1, label = "relaxed")))
  expect_true(out$success)
  expect_identical(out$calibration_table, calibration)
  expect_length(seen, 2L)
  expect_true(all(vapply(seen, identical, logical(1), calibration)))
})

test_that("validate_chronos_output rejects numerically invalid chronos outputs", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  calibration <- data.frame(
    node = as.integer(get_root_num(tree)),
    age.min = 50,
    age.max = 50,
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  ctl <- chronos.control()
  ctl$iter.max <- 3000
  ctl$eval.max <- 3000
  ctl$dual.iter.max <- 50
  out <- suppressWarnings(chronos(tree, lambda = 1, model = "discrete", calibration = calibration, control = ctl, quiet = TRUE))

  expect_true(validate_chronos_output(out)$valid)

  out_bad_ll <- out
  attr(out_bad_ll, "ploglik") <- -1e+100
  expect_false(validate_chronos_output(out_bad_ll)$valid)

  out_zero_edge <- out
  out_zero_edge$edge.length[1] <- 0
  expect_false(validate_chronos_output(out_zero_edge)$valid)
})

test_that("run_chronos_with_restarts treats invalid relaxed chronos result as failure", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "(((t1:0.04671798803,(((t2:1e-08,t3:4.857363743e-08)n6:0.01533183806,t4:0.01287520729)n5:1e-08,t5:0.0001124112221)n4:46.66416385)n3:1.982858307e-05,t6:1e-08)n2:1000,(t7:1e-08,t8:1e-08)n7:1e-08)n1;")
  calibration <- data.frame(
    node = as.integer(c(9, 15, 13, 14, 11, 12)),
    age.min = c(100, 97.47033, 97.47033, 97.47033, 97.47033, 97.47033),
    age.max = c(100, 97.47033, 97.47033, 97.47033, 97.47033, 97.47033),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  ctl <- chronos.control()
  ctl$iter.max <- 1200
  ctl$eval.max <- 1200
  ctl$dual.iter.max <- 80
  out <- suppressWarnings(
    run_chronos_with_restarts(
      phy = tree,
      calibration = calibration,
      chronos_control = ctl,
      chronos_lambda = 1,
      chronos_model = "relaxed",
      context_label = "test-invalid-relaxed",
      max_restarts = 1,
      seed_base = 1
    )
  )
  expect_false(out$success)
  expect_true("try-error" %in% class(out$chronos_out))
  expect_match(as.character(out$chronos_out), "Invalid chronos output", fixed = TRUE)
})

test_that("build_dated_tree_without_chronos returns a bounded ultrametric tree", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  calibration <- data.frame(
    node = as.integer(c(root_num, n1_num)),
    age.min = c(10, 6),
    age.max = c(10, 8),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  out <- build_dated_tree_without_chronos(tree, calibration, root_num)
  expect_true(inherits(out, "phylo"))
  expect_true(is.ultrametric(out))
  ages <- max(node.depth.edgelength(out)) - node.depth.edgelength(out)
  expect_gte(ages[n1_num], 6 - 1e-6)
  expect_lte(ages[n1_num], 8 + 1e-6)
})

test_that("build_dated_tree_without_chronos rejects infeasible lower bounds", {
  tree <- read.tree(text = "((A:1,B:1)n1:1,C:2)root;")
  root_num <- get_root_num(tree)
  n1_num <- get_node_num_by_name(tree, "n1")
  calibration <- data.frame(
    node = as.integer(c(root_num, n1_num)),
    age.min = c(10, 12),
    age.max = c(10, 12),
    soft.bounds = NA,
    stringsAsFactors = FALSE
  )
  out <- build_dated_tree_without_chronos(tree, calibration, root_num)
  expect_s3_class(out, "try-error")
  expect_match(as.character(out), "temporally infeasible")
})
