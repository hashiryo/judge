#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""judge 固有の結果エントリ補強モジュール

Library の build_result_entry が生成する基本エントリに、
judge のダッシュボード・履歴比較で必要なフィールドを付与する。

- time_max_ms:  cases 内の最大実行時間
- memory_max_kb: cases 内の最大メモリ使用量
- cases_hash:   テストケースの同一性判定用ハッシュ

CLI としても使用可能:
  echo '{"cases":[...]}' | python3 enrich_result.py --cases-hash HASH
"""
import hashlib
import json
import sys


def compute_cases_hash(cases: list[dict]) -> str:
    """ケース名からフォールバック用ハッシュを計算する。

    テストケースの .in ファイル内容から計算したハッシュが渡されない場合に使用。
    ケースの追加・削除は検知できるが、内容の変更は検知できない。
    """
    names = sorted(c.get("name", "") for c in cases)
    return hashlib.sha256("\n".join(names).encode()).hexdigest()[:16]


def enrich_result(entry: dict, *, cases_hash: str | None = None) -> dict:
    """結果エントリに judge 固有のフィールドを付与する。"""
    cases = entry.get("cases", [])
    entry["time_max_ms"] = max((c["time_ms"] for c in cases), default=0)
    entry["memory_max_kb"] = max((c["memory_kb"] for c in cases), default=0)
    entry["cases_hash"] = cases_hash if cases_hash else compute_cases_hash(cases)
    return entry


def main() -> None:
    """stdin から JSON を読み、enrich して stdout に書き出す。"""
    cases_hash = None
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--cases-hash" and i + 1 < len(args):
            cases_hash = args[i + 1]
            i += 2
        else:
            i += 1

    entry = json.loads(sys.stdin.read())
    enrich_result(entry, cases_hash=cases_hash)
    sys.stdout.write(json.dumps(entry, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
