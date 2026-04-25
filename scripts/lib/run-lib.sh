#!/usr/bin/env bash
# =============================================================================
# テスト実行共通関数
#
# 使い方: source scripts/lib/run-lib.sh
#
# 呼び出し元が設定すべき変数:
#   ROOT            リポジトリルート (problem.toml 読み取り起点)
#   SCRIPTS_DIR     scripts/ ディレクトリのパス (compare-float-output.py 位置特定用)
#   TC_DIR          testcase キャッシュディレクトリ (get_testcase_dir 用)
#   CUSTOM_TC_DIR   custom testcase ディレクトリ (get_custom_testcase_dir 用)
#
# 任意設定:
#   MLE_MB          メモリ制限 (MB)。未設定なら MLE 判定しない
#   DETAIL_LOG_DIR  エラーログ出力先。未設定ならログ出力しない
#
# 提供関数:
#   parse_problem_toml <problem_dir>   PROBLEM_URL/TLE_SEC/MLE_MB/ERROR_TOL/CALLGRIND_CASE を設定
#   get_testcase_dir <url>             URL の md5 キャッシュディレクトリを echo
#   get_custom_testcase_dir <dir>      custom-testcases/<key> を echo
#   pick_representative_input <hint> <dirs..>  代表テストケース入力 .in を echo
#   append_case_record ...             ケース記録追記 (\x1f 区切り)
#   compile_checker <tc_dir>           checker.cpp のコンパイル
#   run_single_case ...                1 テストケース実行・判定
# =============================================================================

# SCRIPTS_DIR が未設定なら、このファイルの親の親を使う
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# timeout コマンドの解決 (macOS は gtimeout を優先)
if command -v gtimeout >/dev/null 2>&1; then
  TIMEOUT_CMD=gtimeout
elif command -v timeout >/dev/null 2>&1; then
  TIMEOUT_CMD=timeout
else
  echo "Error: 'timeout' (or 'gtimeout') not found. On macOS: brew install coreutils" >&2
  exit 1
fi

# =============================================================================
# problem.toml を解析
# 入力: <problem_dir>  (ROOT からの相対パス)
# 出力変数: PROBLEM_URL / TLE_SEC / MLE_MB / ERROR_TOL / CALLGRIND_CASE
# =============================================================================
parse_problem_toml() {
  local problem_dir="$1"
  local toml_file="${ROOT}/${problem_dir}/problem.toml"
  PROBLEM_URL=""
  TLE_SEC="10"
  MLE_MB="256"
  ERROR_TOL=""
  CALLGRIND_CASE=""

  if [[ ! -f "${toml_file}" ]]; then
    return
  fi

  while IFS= read -r line; do
    if [[ "${line}" =~ ^url\ *=\ *\"([^\"]+)\" ]]; then
      PROBLEM_URL="${BASH_REMATCH[1]}"
    elif [[ "${line}" =~ ^tle\ *=\ *([0-9.]+) ]]; then
      TLE_SEC="${BASH_REMATCH[1]}"
    elif [[ "${line}" =~ ^mle\ *=\ *([0-9]+) ]]; then
      export MLE_MB="${BASH_REMATCH[1]}"
    elif [[ "${line}" =~ ^error\ *=\ *([0-9.eE+-]+) ]]; then
      ERROR_TOL="${BASH_REMATCH[1]}"
    elif [[ "${line}" =~ ^callgrind_case\ *=\ *\"([^\"]+)\" ]]; then
      CALLGRIND_CASE="${BASH_REMATCH[1]}"
    fi
  done < "${toml_file}"
}

# =============================================================================
# testcase ディレクトリ取得 (URL の md5 ベース)
# 入力: <problem_url>
# 出力: 対応する TC_DIR/<md5> を stdout。存在しなければ無出力 + exit 1
# =============================================================================
get_testcase_dir() {
  local problem_url="$1"
  local url_md5
  url_md5=$(echo -n "${problem_url}" | md5sum 2>/dev/null | cut -c1-32 \
    || echo -n "${problem_url}" | md5 -q 2>/dev/null)

  local cache_dir="${TC_DIR}/${url_md5}"
  if [[ -d "${cache_dir}" ]] && [[ "$(ls -A "${cache_dir}" 2>/dev/null)" ]]; then
    echo "${cache_dir}"
    return 0
  fi
  return 1
}

# =============================================================================
# custom testcase ディレクトリ取得 (problem dir を key 化)
# 入力: <problem_dir>
# 出力: CUSTOM_TC_DIR/<key> を stdout。存在しなければ無出力 + exit 1
# =============================================================================
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
# 代表テストケースの入力ファイルを決定
# 入力: <hint_name> <tc_dirs...>
#   hint_name が指定され、かつ <hint>.in が見つかればそれを返す。
#   見つからなければ .in ファイルサイズ最大のものを返す。
# 出力: 入力ファイルパスを stdout
# =============================================================================
pick_representative_input() {
  local hint="$1"
  shift
  local tc_dirs=("$@")

  if [[ -n "${hint}" ]]; then
    for d in "${tc_dirs[@]}"; do
      local candidate="${d}/${hint}.in"
      if [[ -f "${candidate}" ]]; then
        echo "${candidate}"
        return
      fi
    done
    echo "  [WARN] hint=${hint} not found, falling back to largest" >&2
  fi

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

# =============================================================================
# ケース記録の追記 (\x1f 区切り)
# =============================================================================
append_case_record() {
  local records_file="$1"
  local case_name="$2"
  local case_status="$3"
  local case_time="$4"
  local case_mem="$5"
  local case_detail="${6:-}"
  printf '%s\x1f%s\x1f%s\x1f%s\x1f%s\n' \
    "${case_name}" "${case_status}" "${case_time}" "${case_mem}" "${case_detail}" \
    >> "${records_file}"
}

# =============================================================================
# checker のコンパイル
# =============================================================================
compile_checker() {
  local tc_dir="$1"
  local checker_bin="${tc_dir}/checker"

  if [[ ! -f "${tc_dir}/checker.cpp" ]]; then
    echo ""
    return
  fi

  if [[ -x "${checker_bin}" ]]; then
    echo "${checker_bin}"
    return
  fi

  local checker_args=(-std=c++17 -O2)
  [[ -f "${tc_dir}/testlib.h" ]] && checker_args+=("-I${tc_dir}")
  if g++ "${checker_args[@]}" -o "${checker_bin}" "${tc_dir}/checker.cpp" 2>/dev/null; then
    echo "${checker_bin}"
  else
    echo ""
  fi
}

# =============================================================================
# 1つのテストケースを実行して時間・メモリを計測
#
# 引数:
#   $1: binary       実行バイナリ
#   $2: input_file   入力ファイル
#   $3: expected_file 期待出力ファイル
#   $4: tle_sec      タイムアウト秒数
#   $5: error_tolerance 浮動小数点誤差許容値 (空文字なら厳密比較)
#   $6: checker_bin  checker バイナリ (空文字なら不使用)
#
# 環境変数 (任意):
#   case_name       ケース名 (ログファイル名に使用)
#   MLE_MB          メモリ制限 (MB)
#   DETAIL_LOG_DIR  エラーログ出力先
#
# stdout: "STATUS TIME_MS MEMORY_KB [DETAIL]"
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
  local stderr_file
  stderr_file=$(mktemp)

  # 実行 + 計測
  # バイナリの stderr を stderr_file に分離し、/usr/bin/time の出力は time_output へ
  if [[ "$(uname)" == "Darwin" ]]; then
    # shellcheck disable=SC2016  # "$N" は内側の sh -c で展開する
    /usr/bin/time -l sh -c "${TIMEOUT_CMD}"' "$1" "$2" < "$3" > "$4" 2>"$5"' _ \
      "${tle_sec}" "${binary}" "${input_file}" "${output_file}" "${stderr_file}" \
      2>"${time_output}" && true
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
    # shellcheck disable=SC2016  # "$N" は内側の sh -c で展開する
    /usr/bin/time -v sh -c "${TIMEOUT_CMD}"' "$1" "$2" < "$3" > "$4" 2>"$5"' _ \
      "${tle_sec}" "${binary}" "${input_file}" "${output_file}" "${stderr_file}" \
      2>"${time_output}" && true
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

  # 空値のガード
  [[ -z "${elapsed_ms}" ]] && elapsed_ms=0
  [[ -z "${memory_kb}" ]] && memory_kb=0

  # MLE 判定
  if [[ "${status}" == "AC" ]] && [[ -n "${MLE_MB:-}" ]]; then
    local mle_kb=$(( MLE_MB * 1024 ))
    if [[ ${memory_kb} -gt ${mle_kb} ]]; then
      status="MLE"
      detail="used ${memory_kb}KB > limit ${mle_kb}KB"
    fi
  fi

  # 判定 (TLE/RE/MLE でなければ出力を比較)
  if [[ "${status}" == "AC" ]]; then
    if [[ -n "${checker_bin}" ]] && [[ -x "${checker_bin}" ]]; then
      if ! "${checker_bin}" "${input_file}" "${output_file}" "${expected_file}" &>/dev/null; then
        status="WA"
        local actual_head expected_head input_head
        actual_head=$(head -1 "${output_file}" | cut -c1-50)
        expected_head=$(head -1 "${expected_file}" | cut -c1-50)
        input_head=$(head -1 "${input_file}" | cut -c1-80)
        detail="input:[${input_head}] expected:[${expected_head}] actual:[${actual_head}]"
      fi
    elif [[ -n "${error_tolerance}" ]] && [[ "${error_tolerance}" != "0" ]]; then
      if ! python3 "${SCRIPTS_DIR}/compare-float-output.py" \
        --actual "${output_file}" \
        --expected "${expected_file}" \
        --tolerance "${error_tolerance}" >/dev/null 2>&1; then
        status="WA"
        local input_head
        input_head=$(head -1 "${input_file}" | cut -c1-80)
        detail="input:[${input_head}] float compare failed (tolerance=${error_tolerance})"
      fi
    else
      if ! diff <(sed 's/[[:space:]]*$//' "${output_file}") <(sed 's/[[:space:]]*$//' "${expected_file}") &>/dev/null; then
        status="WA"
        local actual_head expected_head input_head
        actual_head=$(head -1 "${output_file}" | cut -c1-50)
        expected_head=$(head -1 "${expected_file}" | cut -c1-50)
        input_head=$(head -1 "${input_file}" | cut -c1-80)
        detail="input:[${input_head}] expected:[${expected_head}] actual:[${actual_head}]"
      fi
    fi
  fi

  # エラー詳細をログディレクトリに保存
  if [[ "${status}" != "AC" ]] && [[ -n "${DETAIL_LOG_DIR:-}" ]]; then
    local log_file="${DETAIL_LOG_DIR}/${case_name:-unknown}.${status}.log"
    {
      echo "Status: ${status}"
      echo "Time: ${elapsed_ms}ms"
      echo "Memory: ${memory_kb}KB"
      if [[ -f "${input_file}" ]]; then
        echo "--- input (first 20 lines) ---"
        head -20 "${input_file}"
      fi
      if [[ "${status}" == "WA" ]] && [[ -f "${output_file}" ]]; then
        echo "--- actual output (first 20 lines) ---"
        head -20 "${output_file}"
        echo "--- expected output (first 20 lines) ---"
        head -20 "${expected_file}"
      fi
      if [[ -s "${stderr_file}" ]]; then
        echo "--- stderr (first 30 lines) ---"
        head -30 "${stderr_file}"
      fi
      echo "--- detail ---"
      echo "${detail}"
    } > "${log_file}" 2>/dev/null
  fi

  rm -f "${output_file}" "${stderr_file}"

  # detail 内の改行や特殊文字をエスケープ
  if [[ -n "${detail}" ]]; then
    detail=$(echo "${detail}" | head -3 | tr '\n' ' ' | sed 's/"/\\"/g' | cut -c1-200)
    echo "${status} ${elapsed_ms} ${memory_kb} ${detail}"
  else
    echo "${status} ${elapsed_ms} ${memory_kb}"
  fi
}
