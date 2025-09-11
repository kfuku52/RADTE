options(repos = c(CRAN="https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
# RADTE が内部で読む/描く際に使うことがある依存
BiocManager::install(c("treeio","ggtree","tidytree"), update=FALSE, ask=FALSE)
install.packages(c("dplyr","phangorn"), quiet=TRUE)
