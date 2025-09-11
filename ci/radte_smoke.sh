#!/usr/bin/env bash
set -euo pipefail

# 走行ログと出力置き場
mkdir -p logs radte_out

echo "[RADTE] version"
radte --version || true

# GeneRax モード（例データ）
echo "[RADTE] GeneRax example"
radte \
  --species_tree data/example_generax_01/species_tree.nwk \
  --generax_nhx data/example_generax_01/gene_tree.nhx \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="$(pwd)/radte_out" \
  --out_prefix=smoke_generax |& tee logs/radte_generax.log

# Notung モード（例データ）
echo "[RADTE] Notung example"
radte \
  --species_tree data/example_notung_01/species_tree.nwk \
  --gene_tree data/example_notung_01/gene_tree.nwk.reconciled \
  --notung_parsable data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="$(pwd)/radte_out" \
  --out_prefix=smoke_notung |& tee logs/radte_notung.log

# ── 最小アサーション（存在・非空・ログ断片） ────────────────────────────
test -s radte_out/smoke_generax_gene_tree_output.nwk
test -s radte_out/smoke_notung_gene_tree_output.nwk

# TSV にヘッダがあること
grep -q $'node\tage' radte_out/*.tsv

# 完走サイン/致命的エラーの簡易チェック
grep -q "Completed: RADTE divergence time estimation" logs/radte_generax.log
grep -q "Completed: RADTE divergence time estimation" logs/radte_notung.log
! grep -qiE "fatal|Execution halted" logs/radte_generax.log
! grep -qiE "fatal|Execution halted" logs/radte_notung.log

echo "[OK] smoke passed"
