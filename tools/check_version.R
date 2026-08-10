#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
against_arg <- args[startsWith(args, "--against=")]
against <- if (length(against_arg) == 0L) NULL else sub("^--against=", "", against_arg[[1]])

extract_version <- function(lines, source) {
  matches <- sub(
    "^radte_version[[:space:]]*=[[:space:]]*'([0-9]+\\.[0-9]+\\.[0-9]+)'[[:space:]]*$",
    "\\1",
    lines[grepl("^radte_version[[:space:]]*=", lines)]
  )
  if (length(matches) != 1L || !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", matches)) {
    stop("Expected exactly one semantic radte_version assignment in ", source)
  }
  matches
}

source_version <- extract_version(readLines(file.path("R", "version.R")), "R/version.R")
bundle_version <- extract_version(readLines("radte.r"), "radte.r")
description <- read.dcf("DESCRIPTION")
description_version <- unname(description[[1, "Version"]])

if (!identical(source_version, bundle_version) || !identical(source_version, description_version)) {
  stop(
    "Version mismatch: R/version.R=", source_version,
    ", radte.r=", bundle_version,
    ", DESCRIPTION=", description_version
  )
}

read_git_file <- function(ref, path) {
  output <- suppressWarnings(system2("git", c("show", shQuote(paste0(ref, ":", path))), stdout = TRUE, stderr = FALSE))
  if (!is.null(attr(output, "status")) && attr(output, "status") != 0L) NULL else output
}

if (!is.null(against)) {
  previous_lines <- read_git_file(against, "R/version.R")
  previous_source <- paste0(against, ":R/version.R")
  if (is.null(previous_lines)) {
    previous_lines <- read_git_file(against, "radte.r")
    previous_source <- paste0(against, ":radte.r")
  }
  if (is.null(previous_lines)) {
    stop("Could not read a version assignment from Git ref: ", against)
  }
  previous_version <- extract_version(previous_lines, previous_source)
  if (package_version(source_version) <= package_version(previous_version)) {
    stop(
      "Version must increase relative to ", against, ": ",
      previous_version, " -> ", source_version
    )
  }
  cat("Version bump verified:", previous_version, "->", source_version, "\n")
} else {
  cat("Version consistency verified:", source_version, "\n")
}
