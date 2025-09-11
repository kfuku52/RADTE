# ci/gen_dummy_trees.R
# needs: ape
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
outdir <- ifelse(length(args) >= 1, args[1], ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## 1) 非ウルトラメトリック種系統樹（敢えて1本だけ短くする）
# 形は ((A,B),(C,D))
non_ultra <- "((A:0.1,B:1):1,(C:1,D:1):1):0;"
writeLines(non_ultra, file.path(outdir, "species_tree.non_ultrametric.nwk"))

## 2) ゼロ長の末端枝を含むがウルトラメトリック
# 形は ((A,B),(C,D)) で root→tip は全て1.0
short_leaves_ultra <- "((A:0,B:0):1,(C:0,D:0):1):0;"
writeLines(short_leaves_ultra, file.path(outdir, "species_tree.short_leaves_ultra.nwk"))

cat("Wrote dummy trees to:", outdir, "\n")
