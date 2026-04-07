#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# ローカルでサンプルテストケースに対してテストを実行する
#
# 使い方:
#   ./scripts/local-test.sh problems/yosupo-unionfind/sol_naive.cpp
#   ./scripts/local-test.sh problems/yosupo-unionfind/   # 全 cpp を実行
# =============================================================================

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# macOS ではデフォルトで clang++ を使う (include/bits/stdc++.h を自前で用意済み)
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"

# ローカルの bits/stdc++.h を include パスに追加
INCLUDE_DIR="${ROOT}/include"
if [[ -d "${INCLUDE_DIR}" ]]; then
  CXXFLAGS="${CXXFLAGS} -I${INCLUDE_DIR}"
fi

if [[ $# -eq 0 ]]; then
  echo "Usage: $0 <cpp-file-or-problem-dir>"
  exit 1
fi

TARGET="$1"

# 対象ファイルを決定
CPP_FILES=()
PROBLEM_DIR=""

if [[ -f "${TARGET}" ]] && [[ "${TARGET}" == *.cpp ]]; then
  CPP_FILES=("${TARGET}")
  PROBLEM_DIR="$(dirname "${TARGET}")"
elif [[ -d "${TARGET}" ]]; then
  PROBLEM_DIR="${TARGET%/}"
  shopt -s nullglob
  for f in "${PROBLEM_DIR}"/*.cpp; do
    CPP_FILES+=("${f}")
  done
  if [[ ${#CPP_FILES[@]} -eq 0 ]]; then
    echo "Error: No .cpp files in ${PROBLEM_DIR}"
    exit 1
  fi
else
  echo "Error: ${TARGET} is not a .cpp file or directory"
  exit 1
fi

# テストケースディレクトリを探す
TC_DIR="${PROBLEM_DIR}/testcases"
if [[ ! -d "${TC_DIR}" ]] || [[ -z "$(ls -A "${TC_DIR}"/*.in 2>/dev/null)" ]]; then
  echo "No testcases found in ${TC_DIR}/"
  echo ""
  echo "Run fetch-samples first:"
  echo "  ./scripts/fetch-samples.sh ${PROBLEM_DIR}/"
  exit 1
fi

# テストケース数を数える
TC_COUNT=$(ls "${TC_DIR}"/*.in 2>/dev/null | wc -l | tr -d ' ')
echo "Testcases: ${TC_DIR}/ (${TC_COUNT} cases)"
echo "Compiler:  ${CXX} ${CXXFLAGS}"
echo ""

# checker の検出・コンパイル
CHECKER_BIN=""
if [[ -f "${TC_DIR}/checker.cpp" ]]; then
  CHECKER_BIN="${TC_DIR}/checker"
  if [[ ! -x "${CHECKER_BIN}" ]] || [[ "${TC_DIR}/checker.cpp" -nt "${CHECKER_BIN}" ]]; then
    echo "Compiling checker..."
    checker_args=(-std=c++17 -O2)
    [[ -f "${TC_DIR}/testlib.h" ]] && checker_args+=("-I${TC_DIR}")
    g++ "${checker_args[@]}" -o "${CHECKER_BIN}" "${TC_DIR}/checker.cpp"
  fi
fi

# 各 cpp ファイルを実行
PASS_ALL=true

for cpp_file in "${CPP_FILES[@]}"; do
  filename=$(basename "${cpp_file}")
  echo "=== ${filename} ==="

  # error tolerance を読み取る
  ERROR_TOL=""
  while IFS= read -r line; do
    if [[ "${line}" =~ //\ *judge-error:\ *([0-9.eE+-]+) ]]; then
      ERROR_TOL="${BASH_REMATCH[1]}"
    fi
  done < <(head -20 "${cpp_file}")

  # コンパイル
  binary=$(mktemp)
  if ! ${CXX} ${CXXFLAGS} -o "${binary}" "${cpp_file}" 2>&1; then
    echo "  CE (Compile Error)"
    echo ""
    PASS_ALL=false
    rm -f "${binary}"
    continue
  fi

  # 各テストケースを実行
  ac_count=0
  total_count=0
  max_time=0

  shopt -s nullglob
  for input_file in "${TC_DIR}"/*.in; do
    [[ -f "${input_file}" ]] || continue
    expected_file="${input_file%.in}.out"
    [[ -f "${expected_file}" ]] || continue

    case_name=$(basename "${input_file}" .in)
    total_count=$((total_count + 1))

    # 実行
    output_file=$(mktemp)
    start_time=$(python3 -c "import time; print(int(time.monotonic_ns()))")

    # timeout コマンド (macOS では gtimeout、なければタイムアウトなし)
    TIMEOUT_CMD=""
    if command -v timeout &>/dev/null; then
      TIMEOUT_CMD="timeout 10"
    elif command -v gtimeout &>/dev/null; then
      TIMEOUT_CMD="gtimeout 10"
    fi

    if ${TIMEOUT_CMD} "${binary}" < "${input_file}" > "${output_file}" 2>/dev/null; then
      end_time=$(python3 -c "import time; print(int(time.monotonic_ns()))")
      elapsed_ms=$(( (end_time - start_time) / 1000000 ))
      [[ ${elapsed_ms} -gt ${max_time} ]] && max_time=${elapsed_ms}

      # 判定
      verdict="AC"
      if [[ -n "${CHECKER_BIN}" ]] && [[ -x "${CHECKER_BIN}" ]]; then
        if ! "${CHECKER_BIN}" "${input_file}" "${output_file}" "${expected_file}" &>/dev/null; then
          verdict="WA"
        fi
      elif [[ -n "${ERROR_TOL}" ]] && [[ "${ERROR_TOL}" != "0" ]]; then
        if ! python3 -c "
import sys
with open('${output_file}') as f: actual = f.read().split()
with open('${expected_file}') as f: expected = f.read().split()
if len(actual) != len(expected): sys.exit(1)
for a, e in zip(actual, expected):
    if abs(float(a) - float(e)) > ${ERROR_TOL}: sys.exit(1)
" 2>/dev/null; then
          verdict="WA"
        fi
      else
        if ! diff <(sed 's/[[:space:]]*$//' "${output_file}") <(sed 's/[[:space:]]*$//' "${expected_file}") &>/dev/null; then
          verdict="WA"
        fi
      fi

      if [[ "${verdict}" == "AC" ]]; then
        echo "  ${case_name}: AC (${elapsed_ms}ms)"
        ac_count=$((ac_count + 1))
      else
        echo "  ${case_name}: WA (${elapsed_ms}ms)"
        echo "    expected: $(head -1 "${expected_file}" | cut -c1-80)"
        echo "    actual:   $(head -1 "${output_file}" | cut -c1-80)"
        PASS_ALL=false
      fi
    else
      exit_code=$?
      end_time=$(python3 -c "import time; print(int(time.monotonic_ns()))")
      elapsed_ms=$(( (end_time - start_time) / 1000000 ))

      if [[ ${exit_code} -eq 124 ]]; then
        echo "  ${case_name}: TLE (>${elapsed_ms}ms)"
      else
        echo "  ${case_name}: RE (exit ${exit_code}, ${elapsed_ms}ms)"
      fi
      PASS_ALL=false
    fi

    rm -f "${output_file}"
  done

  rm -f "${binary}"

  if [[ ${total_count} -eq 0 ]]; then
    echo "  No testcases found"
  else
    echo "  Result: ${ac_count}/${total_count} AC (max ${max_time}ms)"
  fi
  echo ""
done

if ${PASS_ALL}; then
  echo "All passed!"
else
  echo "Some tests failed."
  exit 1
fi
