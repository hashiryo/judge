#!/usr/bin/env python3
"""
変更された問題ディレクトリを検出し、同じ問題の全 .cpp ファイルを含めて出力する。

問題ディレクトリは、problems/ 以下で problem.toml を持つディレクトリとみなす。

出力: JSON (GitHub Actions の matrix 用)
    {"problems": [{"dir": "problems/yosupo-unionfind", "files": ["sol_a.cpp", "sol_b.cpp"]}]}

使い方:
    # push の差分から検出
    python3 scripts/detect-changed.py --before <sha> --after <sha>

    # 全問題を対象
    python3 scripts/detect-changed.py --all
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parent.parent
PROBLEMS_DIR = ROOT / "problems"
EXCLUDED_CPP_DIR_NAMES = {"gen", "testcases"}


def is_submission_cpp(problem_root: Path, cpp_path: Path) -> bool:
    """提出候補に含める .cpp かを判定する"""
    rel = cpp_path.relative_to(problem_root)
    return not any(part in EXCLUDED_CPP_DIR_NAMES for part in rel.parts)


def is_problem_root(path: Path) -> bool:
    """problem.toml を持つディレクトリを問題ルートとみなす"""
    return path.is_dir() and (path / "problem.toml").is_file()


def find_problem_root(path: Path) -> Optional[Path]:
    """指定パスから最も近い問題ルートを親方向に探索する"""
    current = path if path.is_dir() else path.parent
    while True:
        try:
            current.relative_to(PROBLEMS_DIR)
        except ValueError:
            return None

        if is_problem_root(current):
            return current

        if current == PROBLEMS_DIR:
            return None
        current = current.parent


def get_changed_files(before: str, after: str) -> list[str]:
    """git diff で変更されたファイル一覧を取得"""
    try:
        result = subprocess.run(
            ["git", "diff", "--name-only", "--diff-filter=ACMR", f"{before}...{after}"],
            capture_output=True, text=True, check=True, cwd=ROOT,
        )
        return [line.strip() for line in result.stdout.splitlines() if line.strip()]
    except subprocess.CalledProcessError:
        # before が存在しない場合 (初回 push 等) → 全ファイル
        result = subprocess.run(
            ["git", "ls-files"], capture_output=True, text=True, check=True, cwd=ROOT,
        )
        return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def find_problem_dirs_from_files(files: list[str]) -> set[str]:
    """変更ファイルから問題ディレクトリを抽出"""
    dirs = set()
    for f in files:
        file_path = ROOT / f
        try:
            file_path.relative_to(PROBLEMS_DIR)
        except ValueError:
            continue

        problem_root = find_problem_root(file_path)
        if problem_root is not None:
            dirs.add(str(problem_root.relative_to(ROOT)))
    return dirs


def find_all_problem_dirs() -> set[str]:
    """全問題ディレクトリを検出"""
    dirs = set()
    if not PROBLEMS_DIR.exists():
        return dirs
    for toml_file in PROBLEMS_DIR.rglob("problem.toml"):
        child = toml_file.parent
        if child.name.startswith("."):
            continue
        if any(is_submission_cpp(child, cpp) for cpp in child.rglob("*.cpp")):
            dirs.add(str(child.relative_to(ROOT)))
    return dirs


def collect_cpp_files(problem_dir: str) -> list[str]:
    """問題ディレクトリ内の提出候補 .cpp ファイルを収集する"""
    d = ROOT / problem_dir
    files = []
    for cpp in sorted(d.rglob("*.cpp")):
        if is_submission_cpp(d, cpp):
            files.append(str(cpp.relative_to(d)))
    return files


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--before", help="git diff の起点 SHA")
    parser.add_argument("--after", help="git diff の終点 SHA")
    parser.add_argument("--all", action="store_true", help="全問題を対象")
    args = parser.parse_args()

    if args.all:
        problem_dirs = find_all_problem_dirs()
    elif args.before and args.after:
        changed = get_changed_files(args.before, args.after)
        problem_dirs = find_problem_dirs_from_files(changed)
    else:
        print("Error: --all or --before/--after required", file=sys.stderr)
        sys.exit(1)

    problems = []
    for d in sorted(problem_dirs):
        files = collect_cpp_files(d)
        if files:
            problems.append({"dir": d, "files": files})

    output = {"problems": problems}
    matrix_json = json.dumps(output, ensure_ascii=False)

    # GitHub Actions outputs
    github_output = os.environ.get("GITHUB_OUTPUT", "")
    if github_output:
        with open(github_output, "a") as f:
            f.write(f"matrix={matrix_json}\n")
            f.write(f"has-problems={'true' if problems else 'false'}\n")
            # 問題ディレクトリ名のリスト (matrix 用)
            dirs = [p["dir"] for p in problems]
            f.write(f"problem-dirs={json.dumps(dirs)}\n")

    print(matrix_json)
    print(f"Detected {len(problems)} problem(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
