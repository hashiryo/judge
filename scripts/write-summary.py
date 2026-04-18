#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""
テスト結果から GitHub Step Summary 用の比較テーブルを生成する。

- 今回の実行結果を問題×ファイル×環境のテーブルで表示
- history.jsonl から前回の結果を読み込み、差分を表示

出力は stdout に markdown を書き出す。
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULT_DIR = ROOT / ".cache" / "results"
HISTORY_FILE = ROOT / ".results" / "history.jsonl"

ENVS = ["x64-g++", "x64-clang++", "arm-g++", "arm-clang++"]


def load_current_results() -> list[dict]:
    """今回の実行結果を読み込む"""
    results = []
    for rf in sorted(RESULT_DIR.glob("result-*.json")):
        try:
            with open(rf) as f:
                data = json.load(f)
            if isinstance(data, list):
                results.extend(data)
        except Exception:
            pass
    return results


def load_previous_results() -> dict | None:
    """history.jsonl から直前のエントリを読み込む"""
    if not HISTORY_FILE.exists():
        return None
    lines = HISTORY_FILE.read_text().strip().splitlines()
    if len(lines) < 2:
        # 1行しかない = 今回のエントリだけ → 比較対象なし
        return None
    try:
        # 最後から2番目 = 前回
        return json.loads(lines[-2])
    except Exception:
        return None


def format_time(ms: int) -> str:
    if ms >= 1000:
        return f"{ms/1000:.2f}s"
    return f"{ms}ms"


def format_memory(kb: int) -> str:
    if kb >= 1024:
        return f"{kb/1024:.1f}MB"
    return f"{kb}KB"


def format_cell(status: str, time_ms: int, memory_kb: int) -> str:
    if status == "CE":
        return "CE"
    if status == "MLE":
        return f":purple_circle: MLE {format_memory(memory_kb)}"
    icon = {"AC": ":white_check_mark:", "WA": ":x:", "TLE": ":hourglass:", "RE": ":boom:"}.get(status, status)
    return f"{icon} {format_time(time_ms)} / {format_memory(memory_kb)}"


def format_diff(curr_ms: int, prev_ms: int) -> str:
    diff = curr_ms - prev_ms
    if abs(diff) < 3:
        return "="
    sign = "+" if diff > 0 else ""
    return f"{sign}{diff}ms"


def main():
    current = load_current_results()
    if not current:
        print("No results to summarize")
        return

    previous = load_previous_results()

    # 前回結果を (file, env) → result でインデックス化
    prev_index: dict[tuple[str, str], dict] = {}
    if previous:
        for r in previous.get("results", []):
            key = (r.get("file", ""), r.get("environment", ""))
            prev_index[key] = r

    # 問題ごとにグルーピング
    by_problem: dict[str, list[dict]] = defaultdict(list)
    for r in current:
        # 問題ディレクトリを抽出 (ファイルの親ディレクトリ = problem.toml のある場所)
        file_path = r.get("file", "")
        problem_dir = str(Path(file_path).parent)
        if problem_dir.startswith("problems"):
            problem = problem_dir
        else:
            problem = "unknown"
        by_problem[problem].append(r)

    # テーブル生成
    for problem, results in sorted(by_problem.items()):
        print(f"## {problem}\n")

        # ファイル×環境でピボット
        files_data: dict[str, dict[str, dict]] = defaultdict(dict)
        for r in results:
            filename = Path(r["file"]).name
            env = r.get("environment", "")
            files_data[filename][env] = r

        # ヘッダー
        header = "| File |"
        separator = "|:-----|"
        for env in ENVS:
            header += f" {env} |"
            separator += ":---:|"
        print(header)
        print(separator)

        # 各ファイル
        for filename in sorted(files_data.keys()):
            row = f"| `{filename}` |"
            # CE は environment が空の場合がある
            ce_result = files_data[filename].get("")
            for env in ENVS:
                r = files_data[filename].get(env)
                if r is None and ce_result and ce_result.get("status") == "CE":
                    r = ce_result
                if r is None:
                    row += " - |"
                    continue
                cell = format_cell(
                    r.get("status", "?"),
                    r.get("time_max_ms", 0),
                    r.get("memory_max_kb", 0),
                )
                row += f" {cell} |"
            print(row)

        # 前回との差分テーブル
        if prev_index:
            has_diff = False
            for filename in sorted(files_data.keys()):
                for env in ENVS:
                    r = files_data[filename].get(env)
                    if r is None:
                        continue
                    full_path = r.get("file", "")
                    prev = prev_index.get((full_path, env))
                    if prev:
                        has_diff = True
                        break
                if has_diff:
                    break

            if has_diff:
                prev_sha = previous.get("sha_short", "???") if previous else "???"
                print(f"\n### vs previous ({prev_sha})\n")

                header = "| File |"
                separator = "|:-----|"
                for env in ENVS:
                    header += f" {env} |"
                    separator += ":---:|"
                print(header)
                print(separator)

                for filename in sorted(files_data.keys()):
                    row = f"| `{filename}` |"
                    for env in ENVS:
                        r = files_data[filename].get(env)
                        if r is None:
                            row += " - |"
                            continue
                        full_path = r.get("file", "")
                        prev = prev_index.get((full_path, env))
                        if prev is None:
                            row += " (new) |"
                            continue
                        curr_hash = r.get("cases_hash", "")
                        prev_hash = prev.get("cases_hash", "")
                        cases_changed = curr_hash and prev_hash and curr_hash != prev_hash
                        if cases_changed:
                            row += " :recycle: cases changed |"
                        else:
                            curr_ms = r.get("time_max_ms", 0)
                            prev_ms = prev.get("time_max_ms", 0)
                            diff_str = format_diff(curr_ms, prev_ms)
                            curr_status = r.get("status", "?")
                            prev_status = prev.get("status", "?")
                            if curr_status != prev_status:
                                row += f" {prev_status}→{curr_status} ({diff_str}) |"
                            else:
                                row += f" {diff_str} |"
                    print(row)

        print()


if __name__ == "__main__":
    main()
