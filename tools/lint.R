#!/usr/bin/env Rscript

if (!requireNamespace("lintr", quietly = TRUE)) {
  stop("Install the development dependency 'lintr' before running this check.")
}

paths <- c("R", "test", "tools", "benchmark")
lints <- unlist(lapply(paths, function(path) lintr::lint_dir(path)), recursive = FALSE)
if (length(lints) > 0L) {
  print(lints)
  stop("Lint failed with ", length(lints), " finding(s).")
}
cat("Lint passed.\n")
