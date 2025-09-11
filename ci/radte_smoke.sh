#!/usr/bin/env bash
set -euo pipefail

# 作業ディレクトリ
cd "${GITHUB_WORKSPACE}"
WORKDIR="${RUNNER_TEMP:-$PWD}/radte_out"
mkdir -p "$WORKDIR"

echo "[RADTE] versions (R packages)"
Rscript - <<'RS'
cat("R: ", R.Version()$version.string, "\n", sep="")
pkgs <- c("RADTE","ape","treeio","ggtree","tidytree","phangorn")
for (p in pkgs) {
  if (requireNamespace(p, quietly=TRUE)) {
    cat(p, ": ", as.character(utils::packageVersion(p)), "\n", sep="")
  } else {
    cat(p, ": not installed\n", sep="")
  }
}
RS

echo "[RADTE] GeneRax example"
radte \
  --species_tree="${GITHUB_WORKSPACE}/data/example_generax_01/species_tree.nwk" \
  --generax_nhx="${GITHUB_WORKSPACE}/data/example_generax_01/gene_tree.nhx" \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="${WORKDIR}" \
  --out_prefix=smoke_generax

ls -l "${WORKDIR}"/smoke_generax* || true

echo "[RADTE] Notung example"
radte \
  --species_tree="${GITHUB_WORKSPACE}/data/example_notung_01/species_tree.nwk" \
  --gene_tree="${GITHUB_WORKSPACE}/data/example_notung_01/gene_tree.nwk.reconciled" \
  --notung_parsable="${GITHUB_WORKSPACE}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt" \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="${WORKDIR}" \
  --out_prefix=smoke_notung

ls -l "${WORKDIR}"/smoke_notung* || true

echo "[RADTE] done"
