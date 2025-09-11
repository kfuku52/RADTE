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

# 作業パス
WORKDIR="${GITHUB_WORKSPACE}"                     # <- CWD 側にも出るため
ARTDIR="${RUNNER_TEMP:-$PWD}/radte_out"          # <- ここへ集約
mkdir -p "$ARTDIR/generax" "$ARTDIR/notung"

# GeneRax 実行後
# 既知の出力名を片っ端から拾って保全（存在すればコピー）
for f in radte_*.* *_species_tree.* *_gene_tree_output.*; do
  [ -f "$WORKDIR/$f" ] && cp -f "$WORKDIR/$f" "$ARTDIR/generax/" || true
done

# Notung 実行後
for f in radte_*.* *_species_tree.* *_gene_tree_output.*; do
  [ -f "$WORKDIR/$f" ] && cp -f "$WORKDIR/$f" "$ARTDIR/notung/" || true
done

# 仕上げ表示
echo "[ARTIFACTS]"; find "$ARTDIR" -maxdepth 2 -type f -printf '%P\n' || true
