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
project_root <- normalizePath(dirname(radte_path))

# End-to-end examples are relatively expensive and several integration tests
# inspect different artifacts from the same immutable run. Cache successful CLI
# fixtures for the duration of the test process instead of rerunning RADTE for
# every assertion group.
radte_cli_fixture_cache <- new.env(parent = emptyenv())

run_cached_radte_fixture <- function(cache_key, cmd) {
  if (exists(cache_key, envir = radte_cli_fixture_cache, inherits = FALSE)) {
    return(get(cache_key, envir = radte_cli_fixture_cache, inherits = FALSE))
  }
  out_dir <- tempfile(pattern = paste0("radte_cached_", cache_key, "_"))
  dir.create(out_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)
  exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  result <- list(exit_code = exit_code, out_dir = out_dir)
  if (exit_code == 0) {
    assign(cache_key, result, envir = radte_cli_fixture_cache)
  }
  result
}
