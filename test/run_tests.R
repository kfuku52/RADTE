#!/usr/bin/env Rscript

library(testthat)

# Determine the directory of this script
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  # Fallback: assume we are run from project root
  script_dir <- "test"
} else {
  script_dir <- dirname(script_path)
}

test_dir_path <- file.path(script_dir, "testthat")
test_profile <- Sys.getenv("RADTE_TEST_PROFILE", unset = "full")

to_result_df <- function(results) {
  if (is.null(results)) {
    return(data.frame())
  }
  if (is.data.frame(results)) {
    return(results)
  }
  converted <- tryCatch(as.data.frame(results), error = function(e) NULL)
  if (is.null(converted)) {
    return(data.frame())
  }
  converted
}

count_result_col <- function(results, col_name) {
  result_df <- to_result_df(results)
  if (!(col_name %in% names(result_df))) {
    return(0L)
  }
  values <- result_df[[col_name]]
  if (is.logical(values)) {
    return(sum(values, na.rm = TRUE))
  }
  values_num <- suppressWarnings(as.numeric(values))
  return(sum(values_num[is.finite(values_num)], na.rm = TRUE))
}

exit_if_test_failures <- function(results) {
  failed_count <- count_result_col(results, "failed")
  error_count <- count_result_col(results, "error")
  exit_if_fail_count_nonzero(failed_count, error_count)
}

exit_if_fail_count_nonzero <- function(failed_count, error_count) {
  if ((failed_count + error_count) > 0) {
    cat(
      "Test run failed:",
      failed_count, "failed expectation(s),",
      error_count, "error(s).\n"
    )
    quit(save = "no", status = 1)
  }
}

if (test_profile == "full") {
  cat("RADTE test profile: full\n")
  results <- test_dir(
    test_dir_path,
    reporter = "summary"
  )
  exit_if_test_failures(results)
} else if (test_profile == "fast") {
  cat("RADTE test profile: fast\n")
  test_files <- list.files(test_dir_path, pattern = "^test-.*\\.R$", full.names = TRUE)
  fast_exclude <- c(
    "^test-integration\\.R$",
    "^test-integration-extended\\.R$",
    "^test-manual-calculations\\.R$",
    "^test-species-tree-processing\\.R$"
  )
  for (exclude_pattern in fast_exclude) {
    test_files <- test_files[!grepl(exclude_pattern, basename(test_files))]
  }
  if (length(test_files) == 0) {
    stop("No tests selected for fast profile.")
  }
  failed_count <- 0
  error_count <- 0
  for (test_file_path in test_files) {
    file_result <- test_file(
      test_file_path,
      reporter = "summary"
    )
    failed_count <- failed_count + count_result_col(file_result, "failed")
    error_count <- error_count + count_result_col(file_result, "error")
  }
  exit_if_fail_count_nonzero(failed_count, error_count)
} else {
  stop("Unsupported RADTE_TEST_PROFILE: ", test_profile)
}
