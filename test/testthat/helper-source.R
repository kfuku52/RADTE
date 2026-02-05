library(ape)

# Locate radte.r relative to this helper file
# test/testthat/helper-source.R -> ../../radte.r
radte_path <- file.path(dirname(dirname(getwd())), "radte.r")
if (!file.exists(radte_path)) {
  # Fallback: try relative to the test directory
  radte_path <- file.path(getwd(), "..", "..", "radte.r")
}
if (!file.exists(radte_path)) {
  stop("Cannot find radte.r. Expected at: ", radte_path)
}

# Parse radte.r and evaluate only function definitions
# This avoids executing the top-level script logic (library calls, main flow)
exprs <- parse(radte_path)
for (e in exprs) {
  # Match: name = function(...) { ... }
  if (is.call(e) && length(e) >= 3) {
    op <- as.character(e[[1]])
    if (op %in% c("=", "<-")) {
      rhs <- e[[3]]
      if (is.call(rhs) && as.character(rhs[[1]]) == "function") {
        eval(e, envir = globalenv())
      }
    }
  }
}

# Project root for accessing test data
project_root <- normalizePath(file.path(radte_path, ".."))
