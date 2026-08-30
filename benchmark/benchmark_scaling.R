#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  value <- args[startsWith(args, paste0("--", name, "="))]
  if (length(value)) sub(paste0("^--", name, "="), "", value[[1]]) else default
}
root <- option("root", ".")
repetitions <- as.integer(option("repetitions", "5"))
sizes <- args[!startsWith(args, "--")]
sizes <- if (length(sizes)) as.integer(sizes) else c(250L, 500L, 1000L, 2000L, 4000L)
stopifnot(repetitions >= 1L, all(is.finite(sizes)), all(sizes >= 3L))
source(file.path(root, "R", "load.R"))
radte_source_modules(root)

checksum <- function(value) {
  file <- tempfile()
  on.exit(unlink(file))
  saveRDS(value, file, version = 2, compress = FALSE)
  unname(tools::md5sum(file))
}
measure <- function(fun) {
  warmed <- try(fun(), silent = TRUE)
  if (inherits(warmed, "try-error")) return(list(ok = FALSE, median = NA_real_, minimum = NA_real_,
    allocated = NA_real_, hash = NA_character_, error = gsub("[\r\n\t]+", " ", as.character(warmed))))
  times <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    gc()
    times[[i]] <- system.time(value <- fun())[["elapsed"]]
    stopifnot(identical(value, warmed))
  }
  allocated <- NA_real_
  if (capabilities("profmem")) {
    profile <- tempfile()
    Rprofmem(profile)
    tryCatch(fun(), finally = Rprofmem(NULL))
    bytes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(profile))))
    allocated <- sum(bytes, na.rm = TRUE)
    unlink(profile)
  }
  list(ok = TRUE, median = median(times), minimum = min(times), allocated = allocated,
       hash = checksum(warmed), error = "")
}

rows <- list()
for (size in sizes) {
  tree <- ape::stree(size, type = "left")
  tree$tip.label <- paste0("Species", seq_len(size), "_sp")
  tree$node.label <- paste0("node", seq_len(tree$Nnode))
  tree$edge.length <- rep(1, nrow(tree$edge))
  gene <- tree
  gene$tip.label <- paste0(tree$tip.label, "_g1")
  internal <- seq.int(size + 1L, size + tree$Nnode)
  # Strictly feasible constraints, so before/after output equivalence is meaningful.
  ages <- rev(seq_along(internal))
  constraints <- data.frame(gn_node_num = internal, lower_age = ages, upper_age = ages + 0.25)
  root_node <- get_root_num(tree)
  calibration <- data.frame(node = root_node, age.min = size, age.max = size + 1)
  parser <- build_species_parser("legacy")
  workloads <- list(
    clade_match = function() match_internal_nodes_by_clade(tree, tree),
    conflict_scan = function() find_descendant_constraint_conflicts(constraints, tree, root_node),
    constraint_validation = function() stabilize_descendant_constraints(constraints, tree, root_node)$gn_node_table,
    gene_species_mapping = function() map_gene_nodes_to_species_nodes(gene, tree, parser),
    mcmctree_writer = function() build_mcmctree_tree_text(tree, data.frame(), calibration, root_node)
  )
  for (name in names(workloads)) {
    result <- measure(workloads[[name]])
    rows[[length(rows) + 1L]] <- data.frame(tips = size, workload = name, repetitions = repetitions,
      success = result$ok, median_sec = result$median, min_sec = result$minimum,
      allocated_bytes = result$allocated, output_md5 = result$hash, error = result$error)
  }
}
results <- do.call(rbind, rows)
write.table(results, file = option("output", ""), row.names = FALSE, quote = TRUE, sep = "\t")
