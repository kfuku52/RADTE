#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sizes <- if (length(args) == 0L) c(250L, 500L, 1000L, 2000L) else as.integer(args)
source(file.path("R", "load.R"))
radte_source_modules(".")
suppressPackageStartupMessages(library(ape))

results <- do.call(rbind, lapply(sizes, function(size) {
  tree <- stree(size, type = "left")
  internal <- seq.int(Ntip(tree) + 1L, Ntip(tree) + tree$Nnode)
  constraints <- data.frame(
    gn_node_num = internal,
    lower_age = seq_along(internal),
    upper_age = seq_along(internal) + 2
  )
  root <- get_root_num(tree)
  data.frame(
    tips = size,
    clade_match_sec = system.time(match_internal_nodes_by_clade(tree, tree))[["elapsed"]],
    conflict_scan_sec = system.time(
      find_descendant_constraint_conflicts(constraints, tree, root)
    )[["elapsed"]],
    stabilize_sec = system.time(
      stabilize_descendant_constraints(constraints, tree, root)
    )[["elapsed"]]
  )
}))

write.table(results, row.names = FALSE, quote = FALSE, sep = "\t")
