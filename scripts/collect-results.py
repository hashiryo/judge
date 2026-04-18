#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""
各環境のテスト結果 JSON をマージして history.jsonl に1行追記する。

入力: .cache/results/result-*.json (各環境の結果)
出力: .results/history.jsonl に append
"""
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULT_DIR = ROOT / ".cache" / "results"
HISTORY_FILE = ROOT / ".results" / "history.jsonl"


def get_git_info() -> dict:
    """現在の commit 情報を取得"""
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        capture_output=True, text=True, cwd=ROOT,
    ).stdout.strip()

    sha_short = sha[:7] if sha else ""

    ref = os.environ.get("GITHUB_REF", "")
    if not ref:
        ref = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True, text=True, cwd=ROOT,
        ).stdout.strip()

    return {"sha": sha, "sha_short": sha_short, "ref": ref}


def main():
    # 全環境の結果を読み込み
    all_results = []
    result_files = sorted(RESULT_DIR.glob("result-*.json"))

    if not result_files:
        print("No result files found", file=sys.stderr)
        return

    for rf in result_files:
        try:
            with open(rf) as f:
                data = json.load(f)
            if isinstance(data, list):
                all_results.extend(data)
            print(f"  Loaded {rf.name}: {len(data) if isinstance(data, list) else 1} entries")
        except Exception as e:
            print(f"  Error loading {rf.name}: {e}", file=sys.stderr)

    if not all_results:
        print("No results to save", file=sys.stderr)
        return

    # JSONL エントリを構築
    git_info = get_git_info()
    entry = {
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00"),
        "sha": git_info["sha"],
        "sha_short": git_info["sha_short"],
        "ref": git_info["ref"],
        "results": all_results,
    }

    # history.jsonl に追記
    HISTORY_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(HISTORY_FILE, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print(f"Appended to {HISTORY_FILE}: {len(all_results)} results for {git_info['sha_short']}")


if __name__ == "__main__":
    main()
