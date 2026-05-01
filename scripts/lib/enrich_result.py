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
- cpu_model / cpu_arch / compiler_version / cxxflags:
    実行環境情報 (環境変数から取得。GHA ランナー差異を後追いするため)

CLI としても使用可能:
  echo '{"cases":[...]}' | python3 enrich_result.py --cases-hash HASH
"""
import hashlib
import json
import os
import sys


def compute_cases_hash(cases: list[dict]) -> str:
    """ケース名からフォールバック用ハッシュを計算する。

    テストケースの .in ファイル内容から計算したハッシュが渡されない場合に使用。
    ケースの追加・削除は検知できるが、内容の変更は検知できない。
    """
    names = sorted(c.get("name", "") for c in cases)
    return hashlib.sha256("\n".join(names).encode()).hexdigest()[:16]


def _collect_runtime_env() -> dict[str, str]:
    keys = {
        "cpu_model": "RUNNER_CPU_MODEL",
        "cpu_arch": "RUNNER_CPU_ARCH",
        "compiler_version": "COMPILER_VERSION",
        "cxxflags": "CXXFLAGS",
    }
    out: dict[str, str] = {}
    for field, env_name in keys.items():
        value = os.environ.get(env_name, "").strip()
        if value:
            out[field] = value
    return out


def enrich_result(entry: dict, *, cases_hash: str | None = None) -> dict:
    """結果エントリに judge 固有のフィールドを付与する。"""
    cases = entry.get("cases", [])
    entry["time_max_ms"] = max((c["time_ms"] for c in cases), default=0)
    entry["time_total_ms"] = sum((c["time_ms"] for c in cases), start=0)
    entry["memory_max_kb"] = max((c["memory_kb"] for c in cases), default=0)
    entry["cases_hash"] = cases_hash if cases_hash else compute_cases_hash(cases)

    # harness モード由来の純アルゴリズム実行時間 (ns)。
    # 全ケース AC の時のみ集計する。TLE/RE/MLE が混じると ALGO_TIME_NS が
    # 出力されないケースが生じ、max/total が「完走分のみ」になるので
    # 「実態より速く見える」誤解を生む。AC でなければ未付与にして、ダッシュボード
    # 側でも値を出さない (`-` 表示) のが直感に合う。
    if entry.get("status") == "AC":
        algo_times = [c["algo_time_ns"] for c in cases if "algo_time_ns" in c]
        if algo_times:
            entry["algo_time_max_ns"] = max(algo_times)
            entry["algo_time_total_ns"] = sum(algo_times)
    entry.update(_collect_runtime_env())
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

    raw = sys.stdin.read().strip()
    if not raw:
        # 空 stdin (上流の build_result_entry が失敗した等) → 何も書かずに正常終了。
        # JSONDecodeError で CI 全体を巻き込むのを避ける。
        sys.stderr.write("enrich_result: empty stdin, skipping\n")
        return
    entry = json.loads(raw)
    enrich_result(entry, cases_hash=cases_hash)
    sys.stdout.write(json.dumps(entry, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
