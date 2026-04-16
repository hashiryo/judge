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

# Library submodule の共通関数を読み込み
SCRIPTS_DIR="${ROOT}/lib/scripts"
source "${SCRIPTS_DIR}/lib/run-lib.sh"

# macOS ではデフォルトで clang++ を使う (include/bits/stdc++.h を自前で用意済み)
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"

# Library submodule の include を追加 (bits/stdc++.h, debug.hpp)
LIB_INCLUDE="${ROOT}/lib/include"
if [[ -d "${LIB_INCLUDE}" ]]; then
  CXXFLAGS="${CXXFLAGS} -I${LIB_INCLUDE}"
fi

# ARM (macOS等) では SIMDe を有効化
if [[ "$(uname -m)" != "x86_64" ]]; then
  SIMDE_DIR="${LIB_INCLUDE}/simde"
  if [[ -d "${SIMDE_DIR}" ]]; then
    CXXFLAGS="${CXXFLAGS} -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES -I${SIMDE_DIR}"
  fi
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
CHECKER_BIN=$(compile_checker "${TC_DIR}")

# problem.toml から error tolerance を読み取る
ERROR_TOL=""
TOML_FILE="${PROBLEM_DIR}/problem.toml"
if [[ -f "${TOML_FILE}" ]]; then
  while IFS= read -r line; do
    if [[ "${line}" =~ ^error\ *=\ *([0-9.eE+-]+) ]]; then
      ERROR_TOL="${BASH_REMATCH[1]}"
    fi
  done < "${TOML_FILE}"
fi

# 各 cpp ファイルを実行
PASS_ALL=true

for cpp_file in "${CPP_FILES[@]}"; do
  filename=$(basename "${cpp_file}")
  echo "=== ${filename} ==="

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
  max_mem=0

  shopt -s nullglob
  for input_file in "${TC_DIR}"/*.in; do
    [[ -f "${input_file}" ]] || continue
    expected_file="${input_file%.in}.out"
    [[ -f "${expected_file}" ]] || continue

    case_name=$(basename "${input_file}" .in)
    total_count=$((total_count + 1))

    # run_single_case を使用
    result=$(case_name="${case_name}" run_single_case "${binary}" "${input_file}" "${expected_file}" "10" "${ERROR_TOL}" "${CHECKER_BIN}")
    read -r case_status case_time case_mem _ <<< "${result}"

    [[ ${case_time} -gt ${max_time} ]] && max_time=${case_time}
    [[ ${case_mem} -gt ${max_mem} ]] && max_mem=${case_mem}

    if [[ "${case_status}" == "AC" ]]; then
      echo "  ${case_name}: AC (${case_time}ms, ${case_mem}KB)"
      ac_count=$((ac_count + 1))
    else
      echo "  ${case_name}: ${case_status} (${case_time}ms, ${case_mem}KB)"
      if [[ "${case_status}" == "WA" ]]; then
        echo "    expected: $(head -1 "${expected_file}" | cut -c1-80)"
        actual_file=$(mktemp)
        timeout 10 "${binary}" < "${input_file}" > "${actual_file}" 2>/dev/null || true
        echo "    actual:   $(head -1 "${actual_file}" | cut -c1-80)"
        rm -f "${actual_file}"
      fi
      PASS_ALL=false
    fi
  done

  rm -f "${binary}"

  if [[ ${total_count} -eq 0 ]]; then
    echo "  No testcases found"
  else
    echo "  Result: ${ac_count}/${total_count} AC (max ${max_time}ms, ${max_mem}KB)"
  fi
  echo ""
done

if ${PASS_ALL}; then
  echo "All passed!"
else
  echo "Some tests failed."
  exit 1
fi
