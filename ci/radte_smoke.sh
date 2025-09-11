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

# === NEGATIVE / EDGE-CASE TESTS ============================================
echo "[RADTE] generate dummy species trees"
Rscript ci/gen_dummy_trees.R "${RUNNER_TEMP}/radte_dummy"

NON_ULTRA="${RUNNER_TEMP}/radte_dummy/species_tree.non_ultrametric.nwk"
SHORT_ULTRA="${RUNNER_TEMP}/radte_dummy/species_tree.short_leaves_ultra.nwk"

# 例の Notung 入力（付属データ）を使い回し
SP_NOTUNG="${RADTE_WS}/data/example_notung_01/species_tree.nwk"
GN_NOTUNG="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled"
PAR_NOTUNG="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt"

mkdir -p "${RADTE_OUT}/neg_pos"

# --- 1) 非ウルトラメトリック -> 失敗（安全装置で停止）を期待
echo "[RADTE][NEG] non-ultrametric species tree should FAIL with safety stop"
set +e
radte \
  --species_tree="${NON_ULTRA}" \
  --gene_tree="${GN_NOTUNG}" \
  --notung_parsable="${PAR_NOTUNG}" \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --work_dir="${RADTE_OUT}" \
  --out_prefix="neg_non_ultra" \
  > "${RADTE_OUT}/neg_pos/neg_non_ultra.log" 2>&1
rc=$?
set -e

# 結果判定：終了コード≠0 かつ ログに「not ultrametric」＆ stopifnot 由来メッセージが含まれる
if [ $rc -eq 0 ]; then
  echo "[NG] non-ultra tree unexpectedly succeeded"; exit 1
fi
grep -qi "The tree is not ultrametric" "${RADTE_OUT}/neg_pos/neg_non_ultra.log"
grep -qi "sum_adjustment < (sum(tree"     "${RADTE_OUT}/neg_pos/neg_non_ultra.log"
echo "[OK] safety stop for non-ultrametric tree confirmed"

# --- 2) ゼロ長末端枝を含むウルトラメトリック -> パディングされて成功を期待
echo "[RADTE][POS] ultrametric with zero-length terminal branches should PASS with padding"
radte \
  --species_tree="${SHORT_ULTRA}" \
  --gene_tree="${GN_NOTUNG}" \
  --notung_parsable="${PAR_NOTUNG}" \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="${RADTE_OUT}" \
  --out_prefix="pos_short_pad" \
  > "${RADTE_OUT}/neg_pos/pos_short_pad.log" 2>&1

# 成果物の所在は RADTE 側の挙動差に備えて両候補を探す
OUT_NWK=""
for d in "${RADTE_OUT}" "${RADTE_WS}"; do
  f="${d}/radte_gene_tree_output.nwk"
  [ -f "$f" ] && OUT_NWK="$f" && break
done
[ -n "$OUT_NWK" ] || { echo "[NG] output NWK not found"; exit 1; }

# R(ape)で構文・枝長チェック（ログ文言に依存しない）
Rscript - <<RS
suppressPackageStartupMessages(library(ape))
f <- "${OUT_NWK}"
tr <- read.tree(f)
stopifnot(is.list(tr), length(tr$tip.label) > 0)
# 負の枝長がないこと（パディングが効いていれば負値は出ないはず）
stopifnot(all(tr$edge.length >= 0))
cat("OK: parsed", f, "tips=", length(tr$tip.label),
    "min_edge=", min(tr$edge.length), "\n")
RS

echo "[OK] padding flow for short terminal branches confirmed"
