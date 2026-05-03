#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GEN_DIR="${ROOT}/gen"
TESTCASES_DIR="${ROOT}/testcases"

python3 "${GEN_DIR}/make_inputs.py"

REF_BIN="$(mktemp)"
trap 'rm -f "${REF_BIN}"' EXIT

c++ -std=c++17 -O2 -o "${REF_BIN}" "${GEN_DIR}/ref.cpp"

for input_file in "${TESTCASES_DIR}"/*.in; do
  output_file="${input_file%.in}.out"
  "${REF_BIN}" < "${input_file}" > "${output_file}"
done

echo "Generated $(find "${TESTCASES_DIR}" -maxdepth 1 -name '*.in' | wc -l | tr -d ' ') cases"
