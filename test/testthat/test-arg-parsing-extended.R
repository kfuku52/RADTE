# Extended argument parsing tests

test_that("get_parsed_args parses scientific notation as numeric", {
  args <- c("--pad_short_edge=1e-6", "--lambda=1.5e2")
  result <- get_parsed_args(args, print = FALSE)
  expect_true(is.numeric(result$pad_short_edge))
  expect_equal(result$pad_short_edge, 1e-6)
  expect_true(is.numeric(result$lambda))
  expect_equal(result$lambda, 150)
})

test_that("get_parsed_args handles path with special characters", {
  args <- c("--species_tree=/path/to/my tree/file.nwk")
  result <- get_parsed_args(args, print = FALSE)
  expect_equal(result$species_tree, "/path/to/my tree/file.nwk")
})

test_that("get_parsed_args handles NA-like string values", {
  args <- c("--model=discrete", "--value=NA")
  result <- get_parsed_args(args, print = FALSE)
  expect_equal(result$model, "discrete")
  # "NA" stays as string because as.numeric("NA") returns NA, which fails the !is.na check
  expect_equal(result$value, "NA")
})

test_that("get_parsed_args handles many arguments", {
  args <- c("--species_tree=sp.nwk", "--gene_tree=gn.nwk",
            "--notung_parsable=p.txt", "--max_age=1000",
            "--chronos_lambda=1", "--chronos_model=discrete",
            "--pad_short_edge=0.001")
  result <- get_parsed_args(args, print = FALSE)
  expect_equal(length(result), 7)
  expect_equal(result$max_age, 1000)
  expect_equal(result$chronos_lambda, 1)
  expect_equal(result$chronos_model, "discrete")
  expect_equal(result$pad_short_edge, 0.001)
})

test_that("get_parsed_args handles zero value", {
  args <- c("--value=0")
  result <- get_parsed_args(args, print = FALSE)
  expect_true(is.numeric(result$value))
  expect_equal(result$value, 0)
})

test_that("get_parsed_args handles negative numeric value", {
  args <- c("--offset=-10.5")
  result <- get_parsed_args(args, print = FALSE)
  expect_true(is.numeric(result$offset))
  expect_equal(result$offset, -10.5)
})
