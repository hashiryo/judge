#!/usr/bin/env bash
# =============================================================================
# Callgrind 実行スクリプト
#
# 各問題の代表ケース 1 つを Valgrind/Callgrind で実行し、命令数 (Ir) を取得する。
# 結果は .cache/callgrind/callgrind-${ENV_NAME}.jsonl に JSONL 形式で出力。
#
# 代表ケース選定:
#   1. problem.toml の callgrind_case 指定を優先
#   2. なければ testcase ディレクトリ内で .in ファイルサイズ最大のものを使用
#
# 環境変数:
#   CXX, CXXFLAGS, ENV_NAME, PROBLEMS_JSON
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"
ENV_NAME="${ENV_NAME:-local}"
TC_DIR="${ROOT}/.cache/testcases"
CUSTOM_TC_DIR="${ROOT}/.cache/custom-testcases"
OUT_DIR="${ROOT}/.cache/callgrind"
FG_DIR="${ROOT}/.cache/flamegraphs/${ENV_NAME}"
mkdir -p "${OUT_DIR}" "${FG_DIR}"
OUT_JSONL="${OUT_DIR}/callgrind-${ENV_NAME}.jsonl"
: > "${OUT_JSONL}"

# FlameGraph ツールが利用可能か確認 (無ければ SVG 生成はスキップ)
# stackcollapse-callgrind は自作 (scripts/stackcollapse-callgrind.py)、
# flamegraph.pl は Brendan Gregg の FlameGraph リポジトリから取得
FG_TOOLS_DIR="${FG_TOOLS_DIR:-/tmp/FlameGraph}"
STACKCOLLAPSE="${ROOT}/scripts/stackcollapse-callgrind.py"
HAVE_FG=0
if [[ -f "${FG_TOOLS_DIR}/flamegraph.pl" ]] && [[ -f "${STACKCOLLAPSE}" ]]; then
  HAVE_FG=1
fi
echo "FlameGraph: HAVE_FG=${HAVE_FG} (stackcollapse=${STACKCOLLAPSE}, flamegraph=${FG_TOOLS_DIR}/flamegraph.pl)"

if [[ -z "${PROBLEMS_JSON:-}" ]]; then
  echo "Error: PROBLEMS_JSON required"
  exit 1
fi

# testcase ディレクトリ取得 (run-tests.sh と同ロジック)
get_testcase_dir() {
  local url="$1"
  local md5
  md5=$(echo -n "${url}" | md5sum 2>/dev/null | cut -c1-32 || echo -n "${url}" | md5 -q 2>/dev/null)
  local dir="${TC_DIR}/${md5}"
  if [[ -d "${dir}" ]] && [[ "$(ls -A "${dir}" 2>/dev/null)" ]]; then
    echo "${dir}"
  fi
}

get_custom_testcase_dir() {
  local problem_dir="$1"
  local key="${problem_dir//\//_}"
  local dir="${CUSTOM_TC_DIR}/${key}"
  if [[ -d "${dir}" ]] && [[ "$(ls -A "${dir}" 2>/dev/null)" ]]; then
    echo "${dir}"
  fi
}

parse_problem_toml() {
  local toml="${ROOT}/$1/problem.toml"
  PROBLEM_URL=""
  CALLGRIND_CASE=""
  [[ -f "${toml}" ]] || return
  while IFS= read -r line; do
    if [[ "${line}" =~ ^url\ *=\ *\"([^\"]+)\" ]]; then
      PROBLEM_URL="${BASH_REMATCH[1]}"
    elif [[ "${line}" =~ ^callgrind_case\ *=\ *\"([^\"]+)\" ]]; then
      CALLGRIND_CASE="${BASH_REMATCH[1]}"
    fi
  done < "${toml}"
}

# 代表ケースの入力ファイルパスを決定
pick_representative_input() {
  local hint="$1"
  shift
  local tc_dirs=("$@")

  # ヒント指定があれば優先
  if [[ -n "${hint}" ]]; then
    for d in "${tc_dirs[@]}"; do
      local candidate="${d}/${hint}.in"
      if [[ -f "${candidate}" ]]; then
        echo "${candidate}"
        return
      fi
    done
    echo "  [WARN] callgrind_case=${hint} not found, falling back to largest" >&2
  fi

  # サイズ最大の .in を選択
  local largest=""
  local largest_size=0
  for d in "${tc_dirs[@]}"; do
    shopt -s nullglob
    for f in "${d}"/*.in; do
      local sz
      sz=$(stat -c%s "${f}" 2>/dev/null || stat -f%z "${f}" 2>/dev/null || echo 0)
      if [[ ${sz} -gt ${largest_size} ]]; then
        largest_size=${sz}
        largest="${f}"
      fi
    done
  done
  echo "${largest}"
}

echo "Callgrind environment: ${ENV_NAME} (${CXX})"
echo "---"

PROBLEM_COUNT=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; print(len(json.load(sys.stdin)['problems']))")

for i in $(seq 0 $((PROBLEM_COUNT - 1))); do
  PROBLEM_DIR=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; print(json.load(sys.stdin)['problems'][$i]['dir'])")
  FILE_COUNT=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; print(len(json.load(sys.stdin)['problems'][$i]['files']))")

  echo ""
  echo "=== ${PROBLEM_DIR} ==="
  parse_problem_toml "${PROBLEM_DIR}"

  TC_DIRS=()
  if [[ -n "${PROBLEM_URL}" ]]; then
    d=$(get_testcase_dir "${PROBLEM_URL}" || true)
    [[ -n "${d}" ]] && TC_DIRS+=("${d}")
  fi
  d=$(get_custom_testcase_dir "${PROBLEM_DIR}" || true)
  [[ -n "${d}" ]] && TC_DIRS+=("${d}")

  if [[ ${#TC_DIRS[@]} -eq 0 ]]; then
    echo "  [SKIP] no testcases"
    continue
  fi

  INPUT_FILE=$(pick_representative_input "${CALLGRIND_CASE}" "${TC_DIRS[@]}")
  if [[ -z "${INPUT_FILE}" ]] || [[ ! -f "${INPUT_FILE}" ]]; then
    echo "  [SKIP] no representative input"
    continue
  fi
  CASE_NAME=$(basename "${INPUT_FILE}" .in)
  echo "  case: ${CASE_NAME} ($(stat -c%s "${INPUT_FILE}" 2>/dev/null || stat -f%z "${INPUT_FILE}") bytes)"

  for j in $(seq 0 $((FILE_COUNT - 1))); do
    FILENAME=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; print(json.load(sys.stdin)['problems'][$i]['files'][$j])")
    CPP_FILE="${ROOT}/${PROBLEM_DIR}/${FILENAME}"
    REL_PATH="${PROBLEM_DIR}/${FILENAME}"
    [[ -f "${CPP_FILE}" ]] || { echo "  [SKIP] ${REL_PATH}"; continue; }

    BINARY=$(mktemp)
    if ! ${CXX} ${CXXFLAGS} -o "${BINARY}" "${CPP_FILE}" 2>/dev/null; then
      echo "  [CE] ${REL_PATH}"
      rm -f "${BINARY}"
      continue
    fi

    CG_OUT=$(mktemp)
    echo -n "  [CG] ${REL_PATH} ... "
    if valgrind --tool=callgrind --callgrind-out-file="${CG_OUT}" \
        --collect-jumps=no --combine-dumps=yes \
        "${BINARY}" < "${INPUT_FILE}" > /dev/null 2>/dev/null; then
      IR=$(grep -m1 '^totals:' "${CG_OUT}" | awk '{print $2}')
    else
      IR=""
    fi

    FG_REL=""
    if [[ -n "${IR}" ]]; then
      echo -n "${IR}"
      # FlameGraph SVG 生成 (ツールがあれば)
      if [[ ${HAVE_FG} -eq 1 ]]; then
        SVG_PATH="${FG_DIR}/${REL_PATH}.svg"
        mkdir -p "$(dirname "${SVG_PATH}")"
        FOLDED=$(mktemp)
        if python3 "${STACKCOLLAPSE}" "${CG_OUT}" > "${FOLDED}" 2>/dev/null \
            && perl "${FG_TOOLS_DIR}/flamegraph.pl" --title "${REL_PATH} [${ENV_NAME}]" "${FOLDED}" > "${SVG_PATH}" 2>/dev/null; then
          FG_REL="flamegraphs/${ENV_NAME}/${REL_PATH}.svg"
          echo " + FG"
        else
          echo " (FG failed)"
          rm -f "${SVG_PATH}"
        fi
        rm -f "${FOLDED}"
      else
        echo ""
      fi

      python3 -c "
import json, sys
print(json.dumps({
    'file': '${REL_PATH}',
    'environment': '${ENV_NAME}',
    'callgrind_case': '${CASE_NAME}',
    'callgrind_instructions': int('${IR}'),
    'flamegraph_path': '${FG_REL}' or None,
}))" >> "${OUT_JSONL}"
    else
      echo "FAILED"
    fi

    rm -f "${BINARY}" "${CG_OUT}"
  done
done

echo ""
echo "Callgrind results: ${OUT_JSONL} ($(wc -l < "${OUT_JSONL}") entries)"
