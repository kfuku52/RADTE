test_that("get_parsed_args parses key=value pairs correctly", {
  args <- c("--species_tree=/path/to/tree.nwk", "--max_age=100", "--chronos_model=discrete")
  result <- get_parsed_args(args, print = FALSE)
  expect_equal(result$species_tree, "/path/to/tree.nwk")
  expect_equal(result$max_age, 100)
  expect_equal(result$chronos_model, "discrete")
})

test_that("get_parsed_args converts numeric values to numeric type", {
  args <- c("--lambda=1.5", "--integer_val=42")
  result <- get_parsed_args(args, print = FALSE)
  expect_true(is.numeric(result$lambda))
  expect_equal(result$lambda, 1.5)
  expect_true(is.numeric(result$integer_val))
  expect_equal(result$integer_val, 42)
})

test_that("get_parsed_args keeps string values as character", {
  args <- c("--model=discrete", "--path=/some/path")
  result <- get_parsed_args(args, print = FALSE)
  expect_true(is.character(result$model))
  expect_true(is.character(result$path))
})

test_that("get_parsed_args handles single argument", {
  args <- c("--key=value")
  result <- get_parsed_args(args, print = FALSE)
  expect_equal(length(result), 1)
  expect_equal(result$key, "value")
})

test_that("get_parsed_args print=TRUE produces output", {
  args <- c("--key=value")
  expect_output(get_parsed_args(args, print = TRUE), "key")
})
