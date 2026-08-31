#!/usr/bin/env Rscript

# Base R only: this runs before installing any project dependency.
ci_version_only_files <- function(changed, base, head = "HEAD", root = ".") {
  if (!nzchar(base) || grepl("^0+$", base)) return(character())
  candidates <- intersect(changed, c("R/version.R", "radte.r", "DESCRIPTION"))
  Filter(function(path) {
    pattern <- if (path == "DESCRIPTION") {
      "^Version:[[:space:]]+[0-9]+[.][0-9]+[.][0-9]+[[:space:]]*$"
    } else {
      "^radte_version[[:space:]]*=[[:space:]]*'[0-9]+[.][0-9]+[.][0-9]+'[[:space:]]*$"
    }
    status <- system2("git", c("-C", shQuote(root), "diff", "--quiet", "--no-ext-diff", "--no-textconv",
                              paste0("--ignore-matching-lines=", shQuote(pattern)),
                              shQuote(base), shQuote(head), "--", shQuote(path)))
    if (!(status %in% c(0L, 1L))) stop("Could not classify version-only changes in ", path)
    status == 0L
  }, candidates)
}

ci_plan <- function(changed, event, version, version_only = character()) {
  full <- event %in% c("schedule", "workflow_dispatch")
  release <- event == "push" && "R/version.R" %in% changed && grepl("\\.0$", version)
  relevant <- setdiff(changed, version_only)
  code <- full || release || any(grepl(
    "^(R/|test/|tests/|tools/|benchmark/|data/|radte\\.r$|DESCRIPTION$|NAMESPACE$|Makefile$|renv\\.lock$|\\.Rbuildignore$|\\.lintr$|\\.github/)", relevant
  ))
  platforms <- full || release || any(grepl(
    "^(R/|test/|tests/|data/|radte\\.r$|DESCRIPTION$|NAMESPACE$|renv\\.lock$|\\.github/)", relevant
  ))
  matrix <- c(
    '{"os":"ubuntu-latest","r":"4.3","paml":"false"}',
    '{"os":"ubuntu-latest","r":"4.4","paml":"false"}',
    '{"os":"ubuntu-latest","r":"4.5","paml":"false"}',
    '{"os":"ubuntu-latest","r":"4.6","paml":"true"}'
  )
  if (platforms) matrix <- c(matrix,
    '{"os":"macos-latest","r":"4.6","paml":"false"}',
    '{"os":"windows-latest","r":"4.6","paml":"false"}'
  )
  list(code = code, release = release, matrix = paste0('{"include":[', paste(matrix, collapse = ","), "]}"))
}

if (sys.nframe() == 0L) {
  base <- Sys.getenv("RADTE_BASE")
  changed <- if (nzchar(base) && !grepl("^0+$", base)) {
    system2("git", c("diff", "--name-only", shQuote(base), "HEAD", "--"), stdout = TRUE)
  } else "R/version.R"
  if (!is.null(attr(changed, "status"))) stop("Could not determine changed files.")
  version <- read.dcf("DESCRIPTION")[[1, "Version"]]
  version_only <- ci_version_only_files(changed, base)
  plan <- ci_plan(changed, Sys.getenv("GITHUB_EVENT_NAME"), version, version_only)
  lines <- c(paste0("code=", tolower(plan$code)), paste0("release=", tolower(plan$release)),
             paste0("matrix=", plan$matrix))
  cat(paste(lines, collapse = "\n"), "\n", file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
}
