#!/usr/bin/env Rscript

if (!requireNamespace("covr", quietly = TRUE)) {
  stop("Install the development dependency 'covr' before generating coverage.")
}

source_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
test_files <- c(
  file.path("test", "testthat", "helper-source.R"),
  list.files("test/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)
)
test_files <- test_files[!basename(test_files) %in% c(
  "test-integration.R",
  "test-integration-extended.R",
  "test-manual-calculations.R",
  "test-species-tree-processing.R",
  "test-performance-scaling.R"
)]

old_coverage <- Sys.getenv("RADTE_COVERAGE", unset = NA_character_)
Sys.setenv(RADTE_COVERAGE = "true")
on.exit({
  if (is.na(old_coverage)) Sys.unsetenv("RADTE_COVERAGE") else Sys.setenv(RADTE_COVERAGE = old_coverage)
}, add = TRUE)

coverage <- covr::file_coverage(source_files = source_files, test_files = test_files)
summary_text <- capture.output(print(coverage), type = "message")
if (!any(startsWith(summary_text, "Coverage:"))) {
  summary_text <- c(sprintf("Coverage: %.2f%%", covr::percent_coverage(coverage)), summary_text)
}
writeLines(summary_text, "coverage-summary.txt")
covr::to_cobertura(coverage, filename = "coverage.xml")
cat(paste(summary_text, collapse = "\n"), "\n")
