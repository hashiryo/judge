#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# ローカルでサンプルテストケースに対してテストを実行する
#
# 使い方:
#   単一 algo:
#     ./scripts/local-test.sh problems/yosupo-unionfind/algos/naive.hpp
#   ディレクトリ指定 (algos/*.hpp 全て実行):
#     ./scripts/local-test.sh problems/yosupo-unionfind/
# =============================================================================

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# 共通関数を読み込み
SCRIPTS_DIR="${ROOT}/scripts"
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
  echo "Usage: $0 <hpp-file-or-problem-dir>"
  exit 1
fi

TARGET="$1"

# 提出候補 (algos/*.hpp) を ALGOS 配列に集める。
# 各 entry は base.cpp + -DALGO_HPP="algos/xxx.hpp" でコンパイルされる。
ALGOS=()
PROBLEM_DIR=""

if [[ -f "${TARGET}" ]] && [[ "${TARGET}" == *.hpp ]]; then
  algos_dir="$(dirname "${TARGET}")"
  PROBLEM_DIR="$(dirname "${algos_dir}")"
  if [[ "$(basename "${algos_dir}")" != "algos" ]] || [[ ! -f "${PROBLEM_DIR}/base.cpp" ]]; then
    echo "Error: ${TARGET} is not under <problem>/algos/ or base.cpp is missing"
    exit 1
  fi
  ALGOS+=("algos/$(basename "${TARGET}")")
elif [[ -d "${TARGET}" ]]; then
  PROBLEM_DIR="${TARGET%/}"
  if [[ ! -f "${PROBLEM_DIR}/base.cpp" ]] || [[ ! -d "${PROBLEM_DIR}/algos" ]]; then
    echo "Error: ${PROBLEM_DIR} に base.cpp / algos/ がありません"
    exit 1
  fi
  shopt -s nullglob
  for h in "${PROBLEM_DIR}"/algos/*.hpp; do
    name="$(basename "${h}")"
    [[ "${name}" == _* ]] && continue
    ALGOS+=("algos/${name}")
  done
  if [[ ${#ALGOS[@]} -eq 0 ]]; then
    echo "Error: No algos/*.hpp in ${PROBLEM_DIR}"
    exit 1
  fi
else
  echo "Error: ${TARGET} is not a .hpp file or directory"
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

# problem.toml から error tolerance と TLE を読み取る
ERROR_TOL=""
TLE_SEC="10"
TOML_FILE="${PROBLEM_DIR}/problem.toml"
if [[ -f "${TOML_FILE}" ]]; then
  while IFS= read -r line; do
    if [[ "${line}" =~ ^error\ *=\ *([0-9.eE+-]+) ]]; then
      ERROR_TOL="${BASH_REMATCH[1]}"
    elif [[ "${line}" =~ ^tle\ *=\ *([0-9.]+) ]]; then
      TLE_SEC="${BASH_REMATCH[1]}"
    fi
  done < "${TOML_FILE}"
fi

# 各 algo を実行
PASS_ALL=true
SOURCE_CPP="${PROBLEM_DIR}/base.cpp"

for algo_rel in "${ALGOS[@]}"; do
  echo "=== ${algo_rel} ==="

  # コンパイル: base.cpp + -DALGO_HPP="algos/xxx.hpp"
  # base.cpp が "algos/xxx.hpp" を相対 include するため -I${PROBLEM_DIR} を追加
  binary=$(mktemp)
  cmd=("${CXX}")
  # shellcheck disable=SC2206
  cmd+=(${CXXFLAGS})
  cmd+=("-I${PROBLEM_DIR}" "-DALGO_HPP=\"${algo_rel}\"")
  cmd+=(-o "${binary}" "${SOURCE_CPP}")
  if ! "${cmd[@]}" 2>&1; then
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
    result=$(case_name="${case_name}" run_single_case "${binary}" "${input_file}" "${expected_file}" "${TLE_SEC}" "${ERROR_TOL}" "${CHECKER_BIN}")
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
        "${TIMEOUT_CMD}" 10 "${binary}" < "${input_file}" > "${actual_file}" 2>/dev/null || true
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
