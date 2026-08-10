#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
check_only <- identical(args, "--check")
source(file.path("R", "load.R"))
radte_source_modules(".")

target <- file.path("docs", "cli-options.md")
content <- c(
  "# RADTE command-line reference",
  "",
  "This file is generated from `R/cli.R`. Run `make cli-docs` after changing the option schema.",
  "",
  "```text",
  strsplit(render_radte_help(), "\n", fixed = TRUE)[[1]],
  "```"
)

if (check_only) {
  if (!file.exists(target) || !identical(readLines(target, warn = FALSE), content)) {
    stop(target, " is stale. Run: Rscript tools/generate_cli_reference.R")
  }
  cat(target, " is up to date.\n", sep = "")
} else {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  writeLines(content, target, useBytes = TRUE)
  cat("Generated ", target, ".\n", sep = "")
}
