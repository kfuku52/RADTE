#!/usr/bin/env bash
set -euo pipefail

# 追加：ラッパー関数（指定ディレクトリでradteを実行）
run_radte () {
  local outdir="$1"; shift
  mkdir -p "$outdir"
  # CWD を outdir にして実行（--work_dir も合わせて渡す）
  ( cd "$outdir" && radte --work_dir="$outdir" "$@" )
}

# ========== パス定義（成果物はワークスペース直下に統一） ==========
RADTE_WS="${GITHUB_WORKSPACE:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
RUN_TMP="${RUNNER_TEMP:-${RADTE_WS}/_tmp}"

RADTE_OUT="${RADTE_WS}/radte_out"     # ← アーティファクトはここに出す
LOGDIR="${RADTE_WS}/logs"
RADTE_DUMMY="${RUN_TMP}/radte_dummy"

mkdir -p "${RADTE_OUT}/generax" "${RADTE_OUT}/notung" "${LOGDIR}" "${RADTE_DUMMY}"

echo "[RADTE] versions (R packages)"
Rscript - <<'RS'
cat("R: ", R.Version()$version.string, "\n", sep="")
pkgs <- c("RADTE","ape","treeio","ggtree","tidytree","phangorn")
for (p in pkgs) cat(p, ": ", if (requireNamespace(p, quietly=TRUE)) as.character(utils::packageVersion(p)) else "not installed", "\n", sep="")
RS

# すべての実行で共通して使う引数（値は今回のテスト想定）
COMMON_ARGS=(
  --max_age=1000
  --chronos_lambda=1
  --chronos_model=discrete
  --pad_short_edge=0.001
)

# ========== GeneRax ==========
echo "[RADTE] GeneRax example"
run_radte "${RADTE_OUT}/generax" \
  --species_tree="${RADTE_WS}/data/example_generax_01/species_tree.nwk" \
  --generax_nhx="${RADTE_WS}/data/example_generax_01/gene_tree.nhx" \
  "${COMMON_ARGS[@]}" \
  --out_prefix=smoke_generax |& tee "${LOGDIR}/generax.log"

# ========== Notung ==========
echo "[RADTE] Notung example"
run_radte "${RADTE_OUT}/notung" \
  --species_tree="${RADTE_WS}/data/example_notung_01/species_tree.nwk" \
  --gene_tree="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled" \
  --notung_parsable="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt" \
  "${COMMON_ARGS[@]}" \
  --out_prefix=smoke_notung |& tee "${LOGDIR}/notung.log"

# ========== 異常系/境界テスト ==========
echo "[RADTE] generate dummy species trees"
Rscript ci/gen_dummy_trees.R "${RADTE_DUMMY}"

# --- 1) 非ウルトラメトリック → 失敗（終了コード != 0）を期待
echo "[RADTE][NEG] non-ultrametric should FAIL"
mkdir -p "${RADTE_OUT}/neg1"
set +e
( cd "${RADTE_OUT}/neg1" && \
  radte --work_dir="${RADTE_OUT}/neg1" \
    --species_tree="${RADTE_DUMMY}/species_tree.non_ultrametric.nwk" \
    --gene_tree="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled" \
    --notung_parsable="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt" \
    --out_prefix="neg_non_ultra" \
) |& tee "${LOGDIR}/neg_non_ultra.log"
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

# --- 2) ゼロ長枝を含む「遺伝子系統樹」→ パディングされて成功を期待
echo "[RADTE][POS] gene tree with a zero-length edge should PASS (padding)"

# 入力 Notung gene tree をコピーし、最短枝を 0 に書き換えて短枝を強制
MOD_GN="${RUN_TMP}/gene_tree.short.nwk"
GN_IN="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled"
GN_OUT="${MOD_GN}"
mkdir -p "$(dirname "$MOD_GN")"
export GN_IN GN_OUT

echo "GN_IN=${GN_IN}"
ls -l "${GN_IN}"

Rscript - <<'RS'   # ← ここを必ずクオート。bash 変数展開を抑止
suppressPackageStartupMessages(library(ape))
gn <- read.tree(Sys.getenv("GN_IN"))
i  <- which.min(gn$edge.length)
gn$edge.length[i] <- 0  # 最短枝を 0 に
write.tree(gn, Sys.getenv("GN_OUT"))
cat("min_edge_in:", min(gn$edge.length), "\n")
RS

ls -l "${MOD_GN}"

# 実行（種系統樹はオリジナル、遺伝子系統樹だけ改変版を使用）
run_radte "${RADTE_OUT}/pos1" \
  --species_tree="${RADTE_WS}/data/example_notung_01/species_tree.nwk" \
  --gene_tree="${MOD_GN}" \
  --notung_parsable="${RADTE_WS}/data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt" \
  "${COMMON_ARGS[@]}" \
  --out_prefix="pos_short_pad" |& tee "${LOGDIR}/pos1.log"

OUT_NWK="${RADTE_OUT}/pos1/radte_gene_tree_output.nwk"
[ -s "${OUT_NWK}" ] || { 
  echo "❌ no *gene_tree_output.nwk found in ${RADTE_OUT}/pos1"; 
  echo "== pos1.log tail =="
  tail -n 200 "${LOGDIR}/pos1.log" || true
  echo "== search around =="
  find "${RADTE_WS}" -maxdepth 3 -type f -name 'radte_gene_tree_output.nwk' -print || true
  exit 1
}

# 出力ツリーの最短枝が 0.001 以上になっていることを機械検証（パディング確認）
export OUT_NWK
Rscript - <<'RS'    # ← ここもクオート。パスは環境変数で受ける
suppressPackageStartupMessages(library(ape))
f  <- Sys.getenv("OUT_NWK")
tr <- read.tree(f)
stopifnot(is.list(tr), length(tr$tip.label) > 0)
min_e <- min(tr$edge.length)
stopifnot(min_e >= 0.001 - 1e-8)
cat("OK: parsed ", f, ", tips=", length(tr$tip.label),
    ", min_edge=", min_e, " (>= 0.001)\n", sep="")
RS
echo "✅ padding/positive case confirmed"

# ========== 成果物一覧（ログで目視用） ==========
echo "[ARTIFACTS]"
find "${RADTE_OUT}" -maxdepth 2 -type f -printf '%P\n' || true

echo "[RADTE] done"
