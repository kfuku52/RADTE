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

test_dir(
  file.path(script_dir, "testthat"),
  reporter = "summary"
)
