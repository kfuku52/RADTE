#!/usr/bin/env Rscript

# Version/description changes do not invalidate installed dependency caches.
metadata <- read.dcf("DESCRIPTION")
fields <- intersect(c("Depends", "Imports", "Suggests", "LinkingTo", "SystemRequirements"), colnames(metadata))
fingerprint <- c(paste(fields, metadata[1, fields]), readLines("test/install_dependencies.R"))
file <- tempfile()
writeLines(fingerprint, file)
digest <- unname(tools::md5sum(file))
unlink(file)
cat(paste0("dependencies=", digest, "\n"), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)
