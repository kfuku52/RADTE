#!/usr/bin/env bash
set -euo pipefail

# 書き込み可能な作業ディレクトリを用意（RUNNER_TEMP が優先）
WORKDIR="${RUNNER_TEMP:-$GITHUB_WORKSPACE}/radte_out"
mkdir -p "$WORKDIR"

echo "[RADTE] version"
# version / help でも work_dir を必ず渡す
radte --version --work_dir="$WORKDIR"

echo "[RADTE] GeneRax example"
radte \
  --species_tree=data/example_generax_01/species_tree.nwk \
  --generax_nhx=data/example_generax_01/gene_tree.nhx \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="$WORKDIR" \
  --out_prefix=smoke_generax

ls -l "$WORKDIR"/smoke_generax* || true

echo "[RADTE] Notung example"
radte \
  --species_tree=data/example_notung_01/species_tree.nwk \
  --gene_tree=data/example_notung_01/gene_tree.nwk.reconciled \
  --notung_parsable=data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="$WORKDIR" \
  --out_prefix=smoke_notung

ls -l "$WORKDIR"/smoke_notung* || true

echo "[RADTE] done"
