library(ape)
library(testthat)

find_radte_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "R", "load.R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot find the RADTE project root from: ", start)
    }
    current <- parent
  }
}

project_root <- find_radte_project_root()
radte_path <- file.path(project_root, "radte.r")
source(file.path(project_root, "R", "load.R"), local = globalenv())
if (!identical(Sys.getenv("RADTE_COVERAGE"), "true")) {
  radte_source_modules(project_root, envir = globalenv())
}

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
