#!/usr/bin/env Rscript

# Base R only: this runs before installing any project dependency.
ci_plan <- function(changed, event, version) {
  full <- event %in% c("schedule", "workflow_dispatch")
  release <- event == "push" && "R/version.R" %in% changed && grepl("\\.0$", version)
  code <- full || release || any(grepl(
    "^(R/|test/|tests/|tools/|benchmark/|data/|radte\\.r$|DESCRIPTION$|NAMESPACE$|Makefile$|renv\\.lock$|\\.Rbuildignore$|\\.lintr$|\\.github/)", changed
  ))
  platforms <- full || release || any(grepl(
    "^(R/|test/|tests/|data/|radte\\.r$|DESCRIPTION$|NAMESPACE$|renv\\.lock$|\\.github/)", changed
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
  plan <- ci_plan(changed, Sys.getenv("GITHUB_EVENT_NAME"), version)
  lines <- c(paste0("code=", tolower(plan$code)), paste0("release=", tolower(plan$release)),
             paste0("matrix=", plan$matrix))
  cat(paste(lines, collapse = "\n"), "\n", file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
}
