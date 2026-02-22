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

test_that("get_parsed_args returns empty list for empty input", {
  result <- get_parsed_args(character(0), print = FALSE)
  expect_equal(length(result), 0)
})

test_that("get_parsed_args preserves equals signs in argument values", {
  args <- c("--gene_tree=/tmp/a=b.nwk")
  result <- get_parsed_args(args, print = FALSE)
  expect_equal(result$gene_tree, "/tmp/a=b.nwk")
})

test_that("get_parsed_args errors for arguments without equals sign", {
  expect_error(
    get_parsed_args(c("--species_tree"), print = FALSE),
    "key=value"
  )
})

test_that("get_parsed_args errors when argument does not start with --", {
  expect_error(
    get_parsed_args(c("species_tree=sp.nwk"), print = FALSE),
    "key=value"
  )
})

test_that("get_parsed_args errors when key is empty", {
  expect_error(
    get_parsed_args(c("--=sp.nwk"), print = FALSE),
    "key is empty"
  )
})

test_that("get_parsed_args errors when the same key is specified multiple times", {
  expect_error(
    get_parsed_args(c("--max_age=10", "--max_age=20"), print = FALSE),
    "duplicated"
  )
})

test_that("parse_bool_arg accepts common true/false representations", {
  expect_true(parse_bool_arg("true", "--x"))
  expect_true(parse_bool_arg("YES", "--x"))
  expect_false(parse_bool_arg("0", "--x"))
  expect_true(parse_bool_arg(1, "--x"))
  expect_false(parse_bool_arg(FALSE, "--x"))
})

test_that("parse_bool_arg errors on invalid boolean value", {
  expect_error(parse_bool_arg("maybe", "--x"), "boolean")
})

test_that("parse_timeout_arg accepts infinite and zero timeout representations", {
  parsed <- get_parsed_args(c("--timeout=inf"), print = FALSE)
  expect_true(is.infinite(parsed$timeout))
  expect_true(is.infinite(parse_timeout_arg(parsed$timeout, "--timeout")))
  expect_true(is.infinite(parse_timeout_arg("inf", "--timeout")))
  expect_true(is.infinite(parse_timeout_arg("off", "--timeout")))
  expect_true(is.infinite(parse_timeout_arg(Inf, "--timeout")))
  expect_true(is.infinite(parse_timeout_arg(0, "--timeout")))
})

test_that("parse_timeout_arg rejects negative timeout values including -Inf", {
  expect_error(parse_timeout_arg(-1, "--timeout"), "non-negative")
  expect_error(parse_timeout_arg(-Inf, "--timeout"), "non-negative")
})
