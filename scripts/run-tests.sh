#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# テスト実行スクリプト (judge 用)
#
# 問題ディレクトリ内の全 .cpp ファイルをコンパイル・実行し、結果を JSON で出力する。
# 同じ問題の複数提出を同一マシン・同一タイミングで実行して公平に比較する。
#
# 環境変数:
#   CXX        コンパイラ (default: g++)
#   CXXFLAGS   コンパイルフラグ (default: -std=c++17 -O2)
#   ENV_NAME   環境名 (default: local)
#   PROBLEMS_JSON  detect-changed.py の出力 JSON
# =============================================================================

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

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

# =============================================================================
# テストケースディレクトリを取得
# =============================================================================
get_testcase_dir() {
  local problem_url="$1"
  local url_md5
  url_md5=$(echo -n "${problem_url}" | md5sum 2>/dev/null | cut -c1-32 || echo -n "${problem_url}" | md5 -q 2>/dev/null)

  local cache_dir="${TC_DIR}/${url_md5}"
  if [[ -d "${cache_dir}" ]] && [[ "$(ls -A "${cache_dir}" 2>/dev/null)" ]]; then
    echo "${cache_dir}"
    return 0
  fi
  return 1
}

get_custom_testcase_dir() {
  local problem_dir="$1"
  local key="${problem_dir//\//_}"
  local custom_dir="${CUSTOM_TC_DIR}/${key}"
  if [[ -d "${custom_dir}" ]] && [[ "$(ls -A "${custom_dir}" 2>/dev/null)" ]]; then
    echo "${custom_dir}"
    return 0
  fi
  return 1
}

# =============================================================================
# 1つのテストケースを実行して時間・メモリを計測
# =============================================================================
run_single_case() {
  local binary="$1"
  local input_file="$2"
  local expected_file="$3"
  local tle_sec="$4"
  local error_tolerance="$5"
  local checker_bin="${6:-}"
  local output_file
  output_file=$(mktemp)

  local status="AC"
  local elapsed_ms=0
  local memory_kb=0
  local detail=""

  local time_output
  time_output=$(mktemp)

  # 実行 + 計測
  if [[ "$(uname)" == "Darwin" ]]; then
    /usr/bin/time -l timeout "${tle_sec}" "${binary}" < "${input_file}" > "${output_file}" 2>"${time_output}" && true
    local exit_code=$?

    if [[ ${exit_code} -eq 124 ]] || [[ ${exit_code} -eq 137 ]]; then
      status="TLE"
    elif [[ ${exit_code} -ne 0 ]]; then
      status="RE"
      detail="exit code ${exit_code}"
    fi

    elapsed_ms=$(awk '/real/{printf "%.0f", $1 * 1000}' "${time_output}" 2>/dev/null || echo "0")
    memory_kb=$(awk '/maximum resident set size/{printf "%.0f", $1 / 1024}' "${time_output}" 2>/dev/null || echo "0")
  else
    /usr/bin/time -v timeout "${tle_sec}" "${binary}" < "${input_file}" > "${output_file}" 2>"${time_output}" && true
    local exit_code=$?

    if [[ ${exit_code} -eq 124 ]] || [[ ${exit_code} -eq 137 ]]; then
      status="TLE"
    elif [[ ${exit_code} -ne 0 ]]; then
      status="RE"
      detail="exit code ${exit_code}"
    fi

    elapsed_ms=$(grep "Elapsed (wall clock)" "${time_output}" | sed 's/.*: //' | awk -F: '{if (NF==2) printf "%.0f", ($1*60+$2)*1000; else printf "%.0f", $1*1000}' 2>/dev/null || echo "0")
    memory_kb=$(grep "Maximum resident set size" "${time_output}" | awk '{print $NF}' 2>/dev/null || echo "0")
  fi

  rm -f "${time_output}"

  [[ -z "${elapsed_ms}" ]] && elapsed_ms=0
  [[ -z "${memory_kb}" ]] && memory_kb=0

  # 判定
  if [[ "${status}" == "AC" ]]; then
    if [[ -n "${checker_bin}" ]] && [[ -x "${checker_bin}" ]]; then
      if ! "${checker_bin}" "${input_file}" "${output_file}" "${expected_file}" &>/dev/null; then
        status="WA"
      fi
    elif [[ -n "${error_tolerance}" ]] && [[ "${error_tolerance}" != "0" ]]; then
      if ! python3 -c "
import sys
with open('${output_file}') as f: actual = f.read().split()
with open('${expected_file}') as f: expected = f.read().split()
if len(actual) != len(expected): sys.exit(1)
for a, e in zip(actual, expected):
    if abs(float(a) - float(e)) > ${error_tolerance}: sys.exit(1)
" 2>/dev/null; then
        status="WA"
      fi
    else
      if ! diff <(sed 's/[[:space:]]*$//' "${output_file}") <(sed 's/[[:space:]]*$//' "${expected_file}") &>/dev/null; then
        status="WA"
      fi
    fi
  fi

  # エラー詳細をログに保存
  if [[ "${status}" != "AC" ]] && [[ -n "${DETAIL_LOG_DIR:-}" ]]; then
    local log_file="${DETAIL_LOG_DIR}/${case_name:-unknown}.${status}.log"
    {
      echo "Status: ${status}"
      echo "Time: ${elapsed_ms}ms"
      echo "Memory: ${memory_kb}KB"
      if [[ "${status}" == "WA" ]] && [[ -f "${output_file}" ]]; then
        echo "--- actual output (first 20 lines) ---"
        head -20 "${output_file}"
        echo "--- expected output (first 20 lines) ---"
        head -20 "${expected_file}"
      fi
    } > "${log_file}" 2>/dev/null
  fi

  rm -f "${output_file}"

  if [[ -n "${detail}" ]]; then
    detail=$(echo "${detail}" | head -3 | tr '\n' ' ' | sed 's/"/\\"/g' | cut -c1-200)
    echo "${status} ${elapsed_ms} ${memory_kb} ${detail}"
  else
    echo "${status} ${elapsed_ms} ${memory_kb}"
  fi
}

# =============================================================================
# problem.toml を解析
# =============================================================================
parse_problem_toml() {
  local problem_dir="$1"
  local toml_file="${ROOT}/${problem_dir}/problem.toml"
  PROBLEM_URL=""
  TLE_SEC="10"
  ERROR_TOL=""

  if [[ ! -f "${toml_file}" ]]; then
    return
  fi

  while IFS= read -r line; do
    # url = "..."
    if [[ "${line}" =~ ^url\ *=\ *\"([^\"]+)\" ]]; then
      PROBLEM_URL="${BASH_REMATCH[1]}"
    # tle = N
    elif [[ "${line}" =~ ^tle\ *=\ *([0-9.]+) ]]; then
      TLE_SEC="${BASH_REMATCH[1]}"
    # error = N
    elif [[ "${line}" =~ ^error\ *=\ *([0-9.eE+-]+) ]]; then
      ERROR_TOL="${BASH_REMATCH[1]}"
    fi
  done < "${toml_file}"
}

# =============================================================================
# テストケースに対して1つの cpp ファイルを実行
# =============================================================================
run_cpp_file() {
  local cpp_file="$1"
  local tc_dirs=("${@:2}")  # テストケースディレクトリ (複数可)
  local rel_path
  rel_path="${cpp_file#${ROOT}/}"

  # PROBLEM_URL, TLE_SEC, ERROR_TOL は呼び出し元で parse_problem_toml 済み

  # コンパイル
  local binary
  binary=$(mktemp)
  if ! ${CXX} ${CXXFLAGS} -o "${binary}" "${cpp_file}" 2>/dev/null; then
    echo "  [CE] ${rel_path}"
    echo "{\"file\":\"${rel_path}\",\"status\":\"CE\",\"cases\":[]}"
    rm -f "${binary}"
    return
  fi

  echo -n "  [RUN] ${rel_path} "
  local overall_status="AC"
  local cases_json=""
  local case_count=0
  local max_time=0
  local max_mem=0

  for tc_dir in "${tc_dirs[@]}"; do
    # checker の検出・コンパイル
    local checker_bin=""
    if [[ -f "${tc_dir}/checker.cpp" ]]; then
      checker_bin="${tc_dir}/checker"
      if [[ ! -x "${checker_bin}" ]]; then
        local checker_args=(-std=c++17 -O2)
        [[ -f "${tc_dir}/testlib.h" ]] && checker_args+=("-I${tc_dir}")
        g++ "${checker_args[@]}" -o "${checker_bin}" "${tc_dir}/checker.cpp" 2>/dev/null || checker_bin=""
      fi
    fi

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
      local case_status case_time case_mem
      case_status=$(echo "${result}" | awk '{print $1}')
      case_time=$(echo "${result}" | awk '{print $2}')
      case_mem=$(echo "${result}" | awk '{print $3}')

      if [[ "${case_status}" != "AC" ]]; then
        overall_status="${case_status}"
      fi
      [[ ${case_time} -gt ${max_time} ]] && max_time=${case_time}
      [[ ${case_mem} -gt ${max_mem} ]] && max_mem=${case_mem}

      [[ -n "${cases_json}" ]] && cases_json+=","
      cases_json+="{\"name\":\"${case_name}\",\"status\":\"${case_status}\",\"time_ms\":${case_time},\"memory_kb\":${case_mem}}"
      case_count=$((case_count + 1))
    done
  done

  echo "${overall_status} (${case_count} cases, max ${max_time}ms, ${max_mem}KB)"

  echo "{\"file\":\"${rel_path}\",\"problem\":\"${PROBLEM_URL}\",\"environment\":\"${ENV_NAME}\",\"status\":\"${overall_status}\",\"time_max_ms\":${max_time},\"memory_max_kb\":${max_mem},\"cases\":[${cases_json}]}"

  rm -f "${binary}"
}

# =============================================================================
# メイン: 問題ごとに全 cpp を実行
# =============================================================================
EXECUTION_TIME=$(date -u +"%Y-%m-%dT%H:%M:%S+00:00")

echo "Environment: ${ENV_NAME} (${CXX})"
echo "Time: ${EXECUTION_TIME}"
echo "---"

# 結果ファイル (JSON array)
RESULT_FILE="${RESULT_DIR}/result-${ENV_NAME}.json"
echo "[" > "${RESULT_FILE}"
FIRST_ENTRY=true

# PROBLEMS_JSON をパースして問題ごとに処理
PROBLEM_COUNT=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; d=json.load(sys.stdin); print(len(d['problems']))")

for i in $(seq 0 $((PROBLEM_COUNT - 1))); do
  PROBLEM_DIR=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['problems'][$i]['dir'])")
  FILE_COUNT=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; d=json.load(sys.stdin); print(len(d['problems'][$i]['files']))")

  echo ""
  echo "=== ${PROBLEM_DIR} ==="

  # problem.toml から設定を読み込み
  parse_problem_toml "${PROBLEM_DIR}"

  # テストケースディレクトリを収集
  TC_DIRS=()

  if [[ -n "${PROBLEM_URL}" ]]; then
    tc_dir=$(get_testcase_dir "${PROBLEM_URL}" || true)
    if [[ -n "${tc_dir}" ]]; then
      TC_DIRS+=("${tc_dir}")
    fi
  fi

  # 自作テストケース
  custom_dir=$(get_custom_testcase_dir "${PROBLEM_DIR}" || true)
  if [[ -n "${custom_dir}" ]]; then
    TC_DIRS+=("${custom_dir}")
  fi

  if [[ ${#TC_DIRS[@]} -eq 0 ]]; then
    echo "  [SKIP] No testcases available"
    continue
  fi

  # 全 cpp ファイルを実行
  for j in $(seq 0 $((FILE_COUNT - 1))); do
    FILENAME=$(echo "${PROBLEMS_JSON}" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['problems'][$i]['files'][$j])")
    CPP_FILE="${ROOT}/${PROBLEM_DIR}/${FILENAME}"
    if [[ ! -f "${CPP_FILE}" ]]; then
      echo "  [SKIP] ${PROBLEM_DIR}/${FILENAME} (not found)"
      continue
    fi

    result_json=$(run_cpp_file "${CPP_FILE}" "${TC_DIRS[@]}")

    # JSON に追加 (最後の行が JSON)
    json_line=$(echo "${result_json}" | tail -1)
    if [[ "${json_line}" == "{"* ]]; then
      if ${FIRST_ENTRY}; then
        FIRST_ENTRY=false
      else
        echo "," >> "${RESULT_FILE}"
      fi
      echo "${json_line}" >> "${RESULT_FILE}"
    fi
  done
done

echo "]" >> "${RESULT_FILE}"

# 整形
python3 -c "
import json
with open('${RESULT_FILE}') as f:
    data = json.load(f)
with open('${RESULT_FILE}', 'w') as f:
    json.dump(data, f, indent=2)
print(f'Results: ${RESULT_FILE} ({len(data)} entries)')
"
