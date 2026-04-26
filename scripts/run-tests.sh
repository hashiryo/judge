#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# テスト実行スクリプト (judge 用)
#
# 各問題の algos/*.hpp を base.cpp + -DALGO_HPP="..." でビルド・実行し、
# 結果を JSON で出力する。同じ問題の複数提出を同一マシン・同一タイミングで
# 実行して公平に比較する。
#
# 環境変数:
#   CXX        コンパイラ (default: g++)
#   CXXFLAGS   コンパイルフラグ (default: -std=c++17 -O2)
#   ENV_NAME   環境名 (default: local)
#   PROBLEMS_JSON  detect-changed.py の出力 JSON
# =============================================================================

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# 共通関数を読み込み
SCRIPTS_DIR="${ROOT}/scripts"
source "${SCRIPTS_DIR}/lib/run-lib.sh"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"
ENV_NAME="${ENV_NAME:-local}"
TC_DIR="${ROOT}/.cache/testcases"
CUSTOM_TC_DIR="${ROOT}/.cache/custom-testcases"
RESULT_DIR="${ROOT}/.cache/results"

# 引数パース
PROBLEMS_JSON="${PROBLEMS_JSON:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --problems-json) PROBLEMS_JSON="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ -z "${PROBLEMS_JSON}" ]]; then
  echo "Error: --problems-json required or set PROBLEMS_JSON env"
  exit 1
fi

mkdir -p "${RESULT_DIR}"

EXECUTION_TIME=$(date -u +"%Y-%m-%dT%H:%M:%S+00:00")

# parse_problem_toml / get_testcase_dir / get_custom_testcase_dir は run-lib.sh で提供

# =============================================================================
# テストケースに対して 1 つの algos/*.hpp を実行
# =============================================================================
run_cpp_file() {
  local hpp_file="$1"
  local tc_dirs=("${@:2}")
  local rel_path
  rel_path="${hpp_file#${ROOT}/}"

  # base.cpp + -DALGO_HPP="algos/xxx.hpp" でビルド
  local problem_root_abs="${ROOT}/${PROBLEM_DIR}"
  local source_cpp="${problem_root_abs}/base.cpp"
  local algo_rel="${hpp_file#"${problem_root_abs}"/}"
  local -a build_extra=(-DALGO_HPP="\"${algo_rel}\"")

  # コンパイル
  local binary
  binary=$(mktemp)
  local compile_err
  compile_err=$(mktemp)
  if ! ${CXX} ${CXXFLAGS} "${build_extra[@]}" -o "${binary}" "${source_cpp}" 2>"${compile_err}"; then
    echo "  [CE] ${rel_path}" >&2
    head -20 "${compile_err}" | sed 's/^/    /' >&2
    local compile_error_excerpt
    compile_error_excerpt=$(mktemp)
    head -50 "${compile_err}" > "${compile_error_excerpt}"
    python3 "${SCRIPTS_DIR}/collect-run-results.py" build-entry \
      --file "${rel_path}" \
      --problem "${PROBLEM_URL}" \
      --environment "${ENV_NAME}" \
      --status "CE" \
      --last-execution-time "${EXECUTION_TIME}" \
      --compile-error-file "${compile_error_excerpt}" \
      --cases-records /dev/null
    rm -f "${compile_error_excerpt}" "${binary}" "${compile_err}"
    return
  fi
  rm -f "${compile_err}"

  echo -n "  [RUN] ${rel_path} " >&2
  local overall_status="AC"
  local case_count=0
  local case_records_file
  case_records_file=$(mktemp)

  for tc_dir in "${tc_dirs[@]}"; do
    local checker_bin
    checker_bin=$(compile_checker "${tc_dir}")

    shopt -s nullglob
    for input_file in "${tc_dir}"/*.in; do
      [[ -f "${input_file}" ]] || continue

      local expected_file="${input_file%.in}.out"
      [[ -f "${expected_file}" ]] || continue

      local case_name
      case_name=$(basename "${input_file}" .in)

      export DETAIL_LOG_DIR="${ROOT}/.cache/logs/${ENV_NAME}/${rel_path}"
      mkdir -p "${DETAIL_LOG_DIR}" 2>/dev/null || true

      local result
      result=$(case_name="${case_name}" run_single_case "${binary}" "${input_file}" "${expected_file}" "${TLE_SEC}" "${ERROR_TOL}" "${checker_bin}")
      local case_status case_time case_mem case_algo_ns case_detail
      read -r case_status case_time case_mem case_algo_ns case_detail <<< "${result}"

      if [[ "${case_status}" != "AC" ]] && [[ "${overall_status}" == "AC" ]]; then
        overall_status="${case_status}"
      fi

      append_case_record "${case_records_file}" "${case_name}" "${case_status}" "${case_time}" "${case_mem}" "${case_detail:-}" "${case_algo_ns:-0}"
      case_count=$((case_count + 1))
    done
  done

  echo "${overall_status} (${case_count} cases)" >&2

  python3 "${SCRIPTS_DIR}/collect-run-results.py" build-entry \
    --file "${rel_path}" \
    --problem "${PROBLEM_URL}" \
    --environment "${ENV_NAME}" \
    --status "${overall_status}" \
    --last-execution-time "${EXECUTION_TIME}" \
    --cases-records "${case_records_file}"

  rm -f "${case_records_file}" "${binary}"
}

# =============================================================================
# メイン: 問題ごとに全 cpp を実行
# =============================================================================

echo "Environment: ${ENV_NAME} (${CXX})"
echo "Time: ${EXECUTION_TIME}"
echo "---"

# JSONL 中間形式
RESULT_JSONL="${RESULT_DIR}/result-${ENV_NAME}.jsonl"
RESULT_FILE="${RESULT_DIR}/result-${ENV_NAME}.json"
: > "${RESULT_JSONL}"

# PROBLEMS_JSON を TSV 展開: 1 行 = "<dir>\t<file1>\t<file2>..."
while IFS=$'\t' read -r -a PARTS; do
  PROBLEM_DIR="${PARTS[0]}"
  FILES=("${PARTS[@]:1}")

  echo ""
  echo "=== ${PROBLEM_DIR} ==="

  parse_problem_toml "${PROBLEM_DIR}"

  TC_DIRS=()
  if [[ -n "${PROBLEM_URL}" ]]; then
    tc_dir=$(get_testcase_dir "${PROBLEM_URL}" || true)
    [[ -n "${tc_dir}" ]] && TC_DIRS+=("${tc_dir}")
  fi
  custom_dir=$(get_custom_testcase_dir "${PROBLEM_DIR}" || true)
  [[ -n "${custom_dir}" ]] && TC_DIRS+=("${custom_dir}")

  if [[ ${#TC_DIRS[@]} -eq 0 ]]; then
    echo "  [SKIP] No testcases available"
    continue
  fi

  # 問題定義 (テストケース + 計測フレーム + 制約) のハッシュを計算。
  # 内訳: testcases/*.in の内容 + base.cpp (計測ループ) + problem.toml (TLE/MLE 等)。
  # algos/_*.hpp は algo 側の共通ヘッダなので含めない (含めるとリグレッション
  # 比較したい変更で履歴比較がリセットされてしまう)。
  PROBLEM_DEF_FILES=()
  for f in "${ROOT}/${PROBLEM_DIR}/base.cpp" "${ROOT}/${PROBLEM_DIR}/problem.toml"; do
    [[ -f "${f}" ]] && PROBLEM_DEF_FILES+=("${f}")
  done
  CASES_HASH=$(
    {
      find "${TC_DIRS[@]}" -name '*.in' -print0 | sort -z | xargs -0 cat 2>/dev/null
      for f in "${PROBLEM_DEF_FILES[@]}"; do
        printf 'PROBLEM_DEF:%s\n' "${f#"${ROOT}/"}"
        cat "${f}"
      done
    } | shasum -a 256 | cut -c1-16
  )

  for FILENAME in "${FILES[@]}"; do
    CPP_FILE="${ROOT}/${PROBLEM_DIR}/${FILENAME}"
    if [[ ! -f "${CPP_FILE}" ]]; then
      echo "  [SKIP] ${PROBLEM_DIR}/${FILENAME} (not found)"
      continue
    fi

    result_json=$(run_cpp_file "${CPP_FILE}" "${TC_DIRS[@]}")
    if [[ -n "${result_json}" ]]; then
      echo "${result_json}" | python3 "${SCRIPTS_DIR}/lib/enrich_result.py" --cases-hash "${CASES_HASH}" >> "${RESULT_JSONL}"
    fi
  done
done < <(printf '%s' "${PROBLEMS_JSON}" | python3 "${SCRIPTS_DIR}/lib/problems-json.py")

# JSONL → JSON 配列に変換
python3 "${SCRIPTS_DIR}/collect-run-results.py" finalize \
  --in-jsonl "${RESULT_JSONL}" \
  --out-json "${RESULT_FILE}"
