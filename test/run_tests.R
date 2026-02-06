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

if (test_profile == "full") {
  cat("RADTE test profile: full\n")
  test_dir(
    test_dir_path,
    reporter = "summary"
  )
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
  for (test_file_path in test_files) {
    test_file(
      test_file_path,
      reporter = "summary"
    )
  }
} else {
  stop("Unsupported RADTE_TEST_PROFILE: ", test_profile)
}
