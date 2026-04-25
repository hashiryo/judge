#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""PROBLEMS_JSON を TSV に展開する補助スクリプト。

入力: stdin から PROBLEMS_JSON を読み取る。
出力: 1 problem あたり 1 行 TSV、`<dir>\\t<file1>\\t<file2>\\t...` 形式。

bash 側での使い方:
  while IFS=$'\\t' read -r -a PARTS; do
    PROBLEM_DIR="${PARTS[0]}"
    FILES=("${PARTS[@]:1}")
    ...
  done < <(printf '%s' "${PROBLEMS_JSON}" | python3 scripts/lib/problems-json.py)
"""
import json
import sys


def main() -> None:
    data = json.load(sys.stdin)
    for p in data.get("problems", []):
        parts = [p.get("dir", "")] + list(p.get("files", []))
        print("\t".join(parts))


if __name__ == "__main__":
    main()
