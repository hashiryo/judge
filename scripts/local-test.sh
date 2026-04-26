#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# ローカルでサンプルテストケースに対してテストを実行する
#
# 使い方:
#   旧形式:
#     ./scripts/local-test.sh problems/yosupo-unionfind/sol_naive.cpp
#   新形式 (harness):
#     ./scripts/local-test.sh problems/yosupo-unionfind/algos/naive.hpp
#   ディレクトリ指定 (旧/新形式の提出を全て実行):
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
  echo "Usage: $0 <cpp-file-or-problem-dir>"
  exit 1
fi

TARGET="$1"

# 提出候補を決定
# 並列配列:
#   DISPLAY_NAMES[i] : 表示名 (例: "sol_naive.cpp" / "algos/naive.hpp")
#   SOURCES[i]       : コンパイル対象 (.cpp パス。新形式は base.cpp)
#   EXTRA_FLAGS[i]   : 追加フラグ (新形式は -DALGO_HPP=...)
DISPLAY_NAMES=()
SOURCES=()
EXTRA_FLAGS=()
PROBLEM_DIR=""

add_old() {
  local cpp_file="$1"
  DISPLAY_NAMES+=("$(basename "${cpp_file}")")
  SOURCES+=("${cpp_file}")
  EXTRA_FLAGS+=("")
}

add_new() {
  # 新形式: base.cpp を ALGO_HPP=algos/xxx.hpp で再コンパイル
  local prob_dir="$1"
  local hpp_rel="$2"  # 例: algos/naive.hpp
  DISPLAY_NAMES+=("${hpp_rel}")
  SOURCES+=("${prob_dir}/base.cpp")
  # C プリプロセッサに "algos/xxx.hpp" (引用符込み) を渡す必要があるため、
  # シェルレベルで一段エスケープする。
  EXTRA_FLAGS+=("-DALGO_HPP=\"${hpp_rel}\"")
}

collect_dir() {
  local prob_dir="$1"
  shopt -s nullglob
  # 旧形式: <prob_dir>/*.cpp (base.cpp は除外)
  for f in "${prob_dir}"/*.cpp; do
    [[ "$(basename "${f}")" == "base.cpp" ]] && continue
    add_old "${f}"
  done
  # 新形式: base.cpp がある場合のみ algos/*.hpp を列挙 (_ 始まりは除外)
  if [[ -f "${prob_dir}/base.cpp" ]] && [[ -d "${prob_dir}/algos" ]]; then
    for h in "${prob_dir}"/algos/*.hpp; do
      local name
      name="$(basename "${h}")"
      [[ "${name}" == _* ]] && continue
      add_new "${prob_dir}" "algos/${name}"
    done
  fi
}

if [[ -f "${TARGET}" ]] && [[ "${TARGET}" == *.cpp ]]; then
  PROBLEM_DIR="$(dirname "${TARGET}")"
  add_old "${TARGET}"
elif [[ -f "${TARGET}" ]] && [[ "${TARGET}" == *.hpp ]]; then
  # algos/*.hpp 直接指定 (新形式)
  algos_dir="$(dirname "${TARGET}")"
  PROBLEM_DIR="$(dirname "${algos_dir}")"
  if [[ "$(basename "${algos_dir}")" != "algos" ]] || [[ ! -f "${PROBLEM_DIR}/base.cpp" ]]; then
    echo "Error: ${TARGET} is not under <problem>/algos/ or base.cpp is missing"
    exit 1
  fi
  add_new "${PROBLEM_DIR}" "algos/$(basename "${TARGET}")"
elif [[ -d "${TARGET}" ]]; then
  PROBLEM_DIR="${TARGET%/}"
  collect_dir "${PROBLEM_DIR}"
  if [[ ${#SOURCES[@]} -eq 0 ]]; then
    echo "Error: No submission files in ${PROBLEM_DIR}"
    exit 1
  fi
else
  echo "Error: ${TARGET} is not a .cpp/.hpp file or directory"
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

# 各 cpp ファイルを実行
PASS_ALL=true

for i in "${!SOURCES[@]}"; do
  display_name="${DISPLAY_NAMES[$i]}"
  source_file="${SOURCES[$i]}"
  extra="${EXTRA_FLAGS[$i]}"
  echo "=== ${display_name} ==="

  # コンパイル
  # 新形式は base.cpp が "algos/xxx.hpp" を相対参照するため、PROBLEM_DIR を -I に追加
  binary=$(mktemp)
  cmd=("${CXX}")
  # shellcheck disable=SC2206
  cmd+=(${CXXFLAGS})
  cmd+=("-I${PROBLEM_DIR}")
  if [[ -n "${extra}" ]]; then
    cmd+=("${extra}")
  fi
  cmd+=(-o "${binary}" "${source_file}")
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
