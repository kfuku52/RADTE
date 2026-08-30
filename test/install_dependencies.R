#!/usr/bin/env Rscript

cran_repo <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo), timeout = 300)
user_library <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(user_library)) {
  dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_library, .libPaths()))
}

install_cran_with_retry <- function(package, attempts = 3L) {
  for (attempt in seq_len(attempts)) {
    if (requireNamespace(package, quietly = TRUE)) {
      return(invisible(TRUE))
    }
    try(install.packages(package, lib = .libPaths()[[1]]), silent = FALSE)
    if (requireNamespace(package, quietly = TRUE)) {
      return(invisible(TRUE))
    }
  }
  stop("Failed to install ", package, " after ", attempts, " attempts")
}

for (package in c("ape", "testthat", "BiocManager")) {
  install_cran_with_retry(package)
}

install_dev <- tolower(Sys.getenv("RADTE_INSTALL_DEV", unset = "false")) %in%
  c("1", "true", "yes")
if (install_dev) {
  for (package in c("covr", "lintr", "renv")) {
    install_cran_with_retry(package)
  }
}

# R 4.3 selects Bioconductor 3.18 and treeio 1.26.0. That treeio release
# imports tidytree::random_ref, which was removed in tidytree 0.4.7.
if (getRversion() < "4.4.0") {
  install_cran_with_retry("tidytree")
  compatible_tidytree <- "0.4.6"
  archive_url <- paste0(
    sub("/$", "", cran_repo),
    "/src/contrib/Archive/tidytree/tidytree_",
    compatible_tidytree,
    ".tar.gz"
  )
  for (attempt in seq_len(3)) {
    if (packageVersion("tidytree") == compatible_tidytree) {
      break
    }
    try(
      install.packages(archive_url, repos = NULL, type = "source"),
      silent = FALSE
    )
  }
  if (packageVersion("tidytree") != compatible_tidytree) {
    stop("Failed to install compatible tidytree ", compatible_tidytree)
  }
}

for (attempt in seq_len(3)) {
  if (requireNamespace("treeio", quietly = TRUE)) {
    break
  }
  try(
    BiocManager::install("treeio", ask = FALSE, update = FALSE),
    silent = FALSE
  )
}
if (!requireNamespace("treeio", quietly = TRUE)) {
  stop("Failed to install treeio after 3 attempts")
}
