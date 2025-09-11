#!/usr/bin/env bash
set -euo pipefail

# --- 安全な既定パス -----------------------------------------------------------
RADTE_WS="${GITHUB_WORKSPACE:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
RUN_TMP="${RUNNER_TEMP:-${RADTE_WS}/_tmp}"
RADTE_OUT="${RUN_TMP}/radte_out"
RADTE_DUMMY="${RUN_TMP}/radte_dummy"
ARTDIR="${RUN_TMP}/radte_artifacts"
mkdir -p "${RADTE_OUT}" "${RADTE_DUMMY}" "${ARTDIR}/generax" "${ARTDIR}/notung"

echo "[RADTE] versions (R packages)"
Rscript - <<'RS'
cat("R: ", R.Version()$version.string, "\n", sep="")
pkgs <- c("RADTE","ape","treeio","ggtree","tidytree","phangorn")
for (p in pkgs) cat(p, ": ", if (requireNamespace(p, quietly=TRUE)) as.character(utils::packageVersion(p)) else "not installed", "\n", sep="")
RS

# --- GeneRax ------------------------------------------------------------------
echo "[RADTE] GeneRax example"
radte \
  --species_tree="${RADTE_WS}/data/example_generax_01/species_tree.nwk" \
  --generax_nhx="${RADTE_WS}/data/example_generax_01/gene_tree.nhx" \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="${RADTE_OUT}" \
  --out_prefix=smoke_generax

# 成果物コピー（RADTE_OUT から）
shopt -s nullglob
for f in "${RADTE_OUT}"/radte_*.* "${RADTE_OUT}"/*_species_tree.* "${RADTE_OUT}"/*_gene_tree_output.*; do
  cp -f "$f" "${ARTDIR}/generax/"
done
shopt -u nullglob

# --- Notung -------------------------------------------------------------------
echo "[RADTE] Notung example"
radte \
  --species_tree="${RADTE_WS}/data/example_notung_01/species_tree.nwk" \
  --gene_tree="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled" \
  --notung_parsable="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt" \
  --max_age=1000 \
  --chronos_lambda=1 \
  --chronos_model=discrete \
  --pad_short_edge=0.001 \
  --work_dir="${RADTE_OUT}" \
  --out_prefix=smoke_notung

shopt -s nullglob
for f in "${RADTE_OUT}"/radte_*.* "${RADTE_OUT}"/*_species_tree.* "${RADTE_OUT}"/*_gene_tree_output.*; do
  cp -f "$f" "${ARTDIR}/notung/"
done
shopt -u nullglob

# --- 異常系/境界テスト --------------------------------------------------------
echo "[RADTE] generate dummy species trees"
Rscript ci/gen_dummy_trees.R "${RADTE_DUMMY}"

# 1) 非ウルトラメトリック → 失敗を期待（終了コード != 0）
echo "[RADTE][NEG] non-ultrametric should FAIL"
mkdir -p "${RADTE_OUT}/neg1"
set +e
radte \
  --species_tree="${RADTE_DUMMY}/species_tree.non_ultrametric.nwk" \
  --gene_tree="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled" \
  --notung_parsable="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt" \
  --work_dir="${RADTE_OUT}/neg1" \
  --out_prefix="neg_non_ultra" 2>&1 | tee "${RADTE_OUT}/neg1/neg_non_ultra.log"
rc=${PIPESTATUS[0]}
set -e
if [ $rc -eq 0 ]; then
  echo "❌ expected failure (non-ultrametric), but command succeeded"; exit 1
fi
# 失敗時は成果物なしが望ましい
if compgen -G "${RADTE_OUT}/neg1/radte_*" > /dev/null; then
  echo "❌ outputs should not be created on failure"; exit 1
fi
echo "✅ expected failure (non-ultrametric), rc=$rc"

# 2) ゼロ長枝を含む「遺伝子系統樹」→ パディングされて成功を期待
echo "[RADTE][POS] gene tree with a zero-length edge should PASS (padding)"

# 入力 Notung gene tree をコピーし、最短枝を 0 に書き換えて短枝を強制
MOD_GN="${RUN_TMP}/gene_tree.short.nwk"
Rscript - <<RS
suppressPackageStartupMessages(library(ape))
gn <- read.tree("${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled")
i  <- which.min(gn$edge.length)
gn$edge.length[i] <- 0  # 最短枝を 0 に
write.tree(gn, "${MOD_GN}")
cat("min_edge_in:", min(gn$edge.length), "\n")
RS

# 実行（種系統樹はオリジナル、遺伝子系統樹だけ改変版を使用）
mkdir -p "${RADTE_OUT}/pos1"
radte \
  --species_tree="${RADTE_WS}/data/example_notung_01/species_tree.nwk" \
  --gene_tree="${MOD_GN}" \
  --notung_parsable="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt" \
  --pad_short_edge=0.001 \
  --work_dir="${RADTE_OUT}/pos1" \
  --out_prefix="pos_short_pad"

# 代表的成果物をチェック
OUT_NWK="${RADTE_OUT}/pos1/radte_gene_tree_output.nwk"
[ -s "${OUT_NWK}" ] || { echo "❌ expected output not found: ${OUT_NWK}"; exit 1; }

# 出力ツリーの最短枝が 0.001 以上になっていることを機械検証（パディング確認）
Rscript - <<RS
suppressPackageStartupMessages(library(ape))
tr <- read.tree("${OUT_NWK}")
stopifnot(is.list(tr), length(tr$tip.label) > 0)
min_e <- min(tr$edge.length)
stopifnot(min_e >= 0.001 - 1e-8)
cat("OK: parsed ${OUT_NWK}, tips=", length(tr$tip.label),
    ", min_edge=", min_e, " (>= 0.001)\n", sep="")
RS
echo "✅ padding/positive case confirmed"

# --- 成果物一覧（Actions のログで確認用） -------------------------------------
echo "[ARTIFACTS]"
find "${ARTDIR}" -maxdepth 2 -type f -printf '%P\n' || true

echo "[RADTE] done"
