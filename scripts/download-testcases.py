#!/usr/bin/env python3
"""
問題ディレクトリの problem.toml から judge URL を読み取り、テストケースをダウンロードする。
自作テストケース (testcases/) も検出してキャッシュにコピーする。

Library リポジトリの共通ダウンロードモジュール (lib/download.py) を使用。
"""
import json
import os
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TC_CACHE_DIR = ROOT / ".cache" / "testcases"
CUSTOM_TC_DIR = ROOT / ".cache" / "custom-testcases"

YUKICODER_TOKEN = os.environ.get("YUKICODER_TOKEN", "")


def load_download_module():
    """Library submodule の download モジュールを遅延ロードする"""
    scripts_dir = ROOT / "lib" / "scripts"
    sys.path.insert(0, str(scripts_dir))
    from lib.download import DownloadConfig, download_batch

    config = DownloadConfig(
        tc_cache_dir=TC_CACHE_DIR,
        library_checker_dir=ROOT / ".cache" / "library-checker-problems",
        yukicoder_token=YUKICODER_TOKEN,
    )
    return download_batch, config


def parse_problem_toml(problem_dir: Path) -> dict:
    """problem.toml を読み取る (toml パーサーなしの簡易実装)"""
    toml_path = problem_dir / "problem.toml"
    config: dict[str, str] = {}
    if not toml_path.exists():
        return config
    for line in toml_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = re.match(r'(\w+)\s*=\s*"([^"]*)"', line)
        if m:
            config[m.group(1)] = m.group(2)
            continue
        m = re.match(r"(\w+)\s*=\s*(\S+)", line)
        if m:
            config[m.group(1)] = m.group(2)
    return config


def collect_problem_urls(problems_json: dict) -> dict[str, list[str]]:
    """問題 JSON から URL → ファイルリストのマッピングを作成"""
    urls: dict[str, list[str]] = {}
    for problem in problems_json.get("problems", []):
        problem_dir = ROOT / problem["dir"]
        config = parse_problem_toml(problem_dir)
        url = config.get("url")
        if url:
            if url not in urls:
                urls[url] = []
            for filename in problem["files"]:
                urls[url].append(
                    str((problem_dir / filename).relative_to(ROOT))
                )
    return urls


def collect_custom_testcases(
    problems_json: dict,
) -> list[tuple[str, Path]]:
    """自作テストケースを持つ問題ディレクトリを検出"""
    result: list[tuple[str, Path]] = []
    for problem in problems_json.get("problems", []):
        problem_dir = ROOT / problem["dir"]
        tc_dir = problem_dir / "testcases"
        if tc_dir.is_dir():
            in_files = list(tc_dir.glob("*.in"))
            if in_files:
                result.append((problem["dir"], tc_dir))
    return result


def setup_custom_testcases(problem_dir_name: str, tc_dir: Path) -> None:
    """自作テストケースをキャッシュディレクトリにコピー"""
    cache_dir = CUSTOM_TC_DIR / problem_dir_name.replace("/", "_")
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for in_file in sorted(tc_dir.glob("*.in")):
        out_file = tc_dir / f"{in_file.stem}.out"
        if out_file.exists():
            shutil.copy2(in_file, cache_dir / in_file.name)
            shutil.copy2(out_file, cache_dir / out_file.name)
            count += 1

    checker = tc_dir / "checker.cpp"
    if checker.exists():
        shutil.copy2(checker, cache_dir / "checker.cpp")

    print(f"  Custom testcases for {problem_dir_name}: {count} cases")


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--problems-json",
        required=True,
        help="detect-changed.py の出力 JSON",
    )
    args = parser.parse_args()

    problems_json = json.loads(args.problems_json)
    urls = collect_problem_urls(problems_json)
    print(f"Found {len(urls)} unique problem URLs")

    # 自作テストケースの処理 (judge 固有)
    custom_tcs = collect_custom_testcases(problems_json)
    for problem_dir_name, tc_dir in custom_tcs:
        setup_custom_testcases(problem_dir_name, tc_dir)

    if not urls:
        if custom_tcs:
            print("Only custom testcases, no downloads needed")
        else:
            print("No problem URLs found")
        return

    # 共通ダウンロードモジュールを使用
    download_batch, config = load_download_module()
    summary = download_batch(urls, config)

    print(f"\nTestcase download summary:")
    print(f"  Cached:     {summary.cached}")
    print(f"  Downloaded: {summary.downloaded}")
    print(f"  Failed:     {summary.failed}")
    print(f"  Total URLs: {summary.total}")
    print(f"  Custom:     {len(custom_tcs)} problem(s)")


if __name__ == "__main__":
    main()
