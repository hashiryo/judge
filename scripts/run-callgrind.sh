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
SCRIPTS_DIR="${ROOT}/scripts"
# shellcheck source=lib/run-lib.sh
source "${SCRIPTS_DIR}/lib/run-lib.sh"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"
ENV_NAME="${ENV_NAME:-local}"
TC_DIR="${ROOT}/.cache/testcases"
CUSTOM_TC_DIR="${ROOT}/.cache/custom-testcases"
OUT_DIR="${ROOT}/.cache/callgrind"
mkdir -p "${OUT_DIR}"
OUT_JSONL="${OUT_DIR}/callgrind-${ENV_NAME}.jsonl"
: > "${OUT_JSONL}"

if [[ -z "${PROBLEMS_JSON:-}" ]]; then
  echo "Error: PROBLEMS_JSON required"
  exit 1
fi

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
        --cache-sim=yes --branch-sim=yes \
        --dump-instr=yes --compress-pos=no \
        "${BINARY}" < "${INPUT_FILE}" > /dev/null 2>/dev/null; then
      EVENTS_LINE=$(grep -m1 '^events:' "${CG_OUT}" | sed 's/^events: //')
      TOTALS_LINE=$(grep -m1 '^totals:' "${CG_OUT}" | sed 's/^totals: //')
    else
      EVENTS_LINE=""
      TOTALS_LINE=""
    fi

    if [[ -n "${TOTALS_LINE}" ]]; then
      # SIMD 命令の Ir 集計 (x86: v prefix, arm: vec reg)
      # 出力: "<simd_ir> <binary_ir>" (binary_ir はライブラリ除くバイナリ本体の Ir)
      SIMD_COUNTS=$(python3 "${ROOT}/scripts/count-simd.py" "${BINARY}" "${CG_OUT}" 2>/dev/null || echo "0 0")
      SIMD_IR=$(echo "${SIMD_COUNTS}" | awk '{print $1}')
      BINARY_IR=$(echo "${SIMD_COUNTS}" | awk '{print $2}')

      JSON_LINE=$(EVENTS="${EVENTS_LINE}" TOTALS="${TOTALS_LINE}" \
        SIMD_IR="${SIMD_IR}" BINARY_IR="${BINARY_IR}" \
        REL_PATH="${REL_PATH}" ENV_NAME="${ENV_NAME}" CASE_NAME="${CASE_NAME}" \
        python3 -c '
import json, os
events = os.environ["EVENTS"].split()
totals = list(map(int, os.environ["TOTALS"].split()))
m = dict(zip(events, totals))
simd_ir = int(os.environ["SIMD_IR"])
binary_ir = int(os.environ["BINARY_IR"])
# SIMD ratio: バイナリ本体 Ir に対する SIMD 命令の割合 (ライブラリコード除外)
simd_ratio = (simd_ir / binary_ir) if binary_ir > 0 else 0.0
print(json.dumps({
    "file": os.environ["REL_PATH"],
    "environment": os.environ["ENV_NAME"],
    "callgrind_case": os.environ["CASE_NAME"],
    "callgrind_instructions": m.get("Ir", 0),
    "callgrind_d1_misses": m.get("D1mr", 0) + m.get("D1mw", 0),
    "callgrind_ll_misses": m.get("DLmr", 0) + m.get("DLmw", 0),
    "callgrind_branch_misses": m.get("Bcm", 0),
    "callgrind_simd_instructions": simd_ir,
    "callgrind_binary_instructions": binary_ir,
    "callgrind_simd_ratio": simd_ratio,
}))')
      echo "${JSON_LINE}" >> "${OUT_JSONL}"
      IR=$(echo "${TOTALS_LINE}" | awk '{print $1}')
      echo "Ir=${IR} SIMD=${SIMD_IR}/${BINARY_IR}"
    else
      echo "FAILED"
    fi

    rm -f "${BINARY}" "${CG_OUT}"
  done
done

echo ""
echo "Callgrind results: ${OUT_JSONL} ($(wc -l < "${OUT_JSONL}") entries)"
