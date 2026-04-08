#!/usr/bin/env python3
"""
問題のサンプルケースをダウンロードして testcases/ に保存する。

使い方:
  python3 scripts/fetch-samples.py problems/yosupo-unionfind/
  python3 scripts/fetch-samples.py problems/yosupo-unionfind/sol_naive.cpp

対応サービス:
  - yosupo judge: GitHub の example_*.in + sol/correct.cpp で出力生成
  - AOJ: サンプル専用 API (judgedat.u-aizu.ac.jp/testcases/samples/)
  - yukicoder: 問題ページ HTML からスクレイピング
"""
import html as html_mod
import json
import re
import subprocess
import sys
import tempfile
import urllib.request
import urllib.error
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def parse_problem_toml(problem_dir: Path) -> dict:
    """problem.toml を読み取る"""
    toml_path = problem_dir / "problem.toml"
    config = {}
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
        m = re.match(r'(\w+)\s*=\s*(\S+)', line)
        if m:
            config[m.group(1)] = m.group(2)
    return config


# ============================================================
# yosupo judge (Library Checker)
# ============================================================

def fetch_yosupo(url: str, out_dir: Path) -> int:
    """
    yosupo judge のサンプルを取得。
    library-checker-problems の example_*.in と sol/correct.cpp を使う。
    """
    m = re.search(r"/problem/(\w+)", url)
    if not m:
        print(f"  Cannot parse problem ID from {url}", file=sys.stderr)
        return 0
    problem_id = m.group(1)

    # GitHub API でリポジトリのツリーから問題ディレクトリのパスを特定
    print(f"  Finding {problem_id} in library-checker-problems...")
    tree_url = "https://api.github.com/repos/yosupo06/library-checker-problems/git/trees/master?recursive=1"
    try:
        req = urllib.request.Request(tree_url)
        req.add_header("User-Agent", "judge-fetch-samples")
        with urllib.request.urlopen(req, timeout=30) as resp:
            tree = json.loads(resp.read())
    except Exception as e:
        print(f"  GitHub API error: {e}", file=sys.stderr)
        return 0

    # 問題ディレクトリと example ファイルを探す
    problem_prefix = None
    example_files = []
    sol_correct_path = None

    for item in tree.get("tree", []):
        path = item["path"]
        if path.endswith(f"/{problem_id}/info.toml"):
            problem_prefix = path.rsplit("/info.toml", 1)[0]

    if not problem_prefix:
        print(f"  Problem {problem_id} not found in repo", file=sys.stderr)
        return 0

    for item in tree.get("tree", []):
        path = item["path"]
        if path.startswith(problem_prefix + "/"):
            rel = path[len(problem_prefix) + 1:]
            if rel.startswith("gen/example_") and rel.endswith(".in"):
                example_files.append(path)
            elif rel == "sol/correct.cpp":
                sol_correct_path = path

    if not example_files:
        print(f"  No example files found for {problem_id}", file=sys.stderr)
        return 0

    if not sol_correct_path:
        print(f"  No sol/correct.cpp found for {problem_id}", file=sys.stderr)
        return 0

    print(f"  Found {len(example_files)} example(s), downloading...")

    base_raw = "https://raw.githubusercontent.com/yosupo06/library-checker-problems/master"
    in_files = []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        for ex_path in sorted(example_files):
            filename = Path(ex_path).name
            raw_url = f"{base_raw}/{ex_path}"
            try:
                with urllib.request.urlopen(raw_url, timeout=30) as resp:
                    data = resp.read().decode()
                (tmpdir / filename).write_text(data)
                in_files.append(filename)
            except Exception as e:
                print(f"  Failed to download {filename}: {e}", file=sys.stderr)

        if not in_files:
            return 0

        # sol/correct.cpp をダウンロードしてコンパイル
        print(f"  Compiling sol/correct.cpp...")
        sol_url = f"{base_raw}/{sol_correct_path}"
        try:
            with urllib.request.urlopen(sol_url, timeout=30) as resp:
                sol_code = resp.read().decode()
            (tmpdir / "correct.cpp").write_text(sol_code)
        except Exception as e:
            print(f"  Failed to download correct.cpp: {e}", file=sys.stderr)
            return 0

        # include/bits/stdc++.h があれば -I に追加
        include_dir = ROOT / "include"
        compile_flags = ["-std=c++17", "-O2"]
        if include_dir.is_dir():
            compile_flags.append(f"-I{include_dir}")

        binary = tmpdir / "correct"
        for cxx in ["g++", "c++"]:
            try:
                subprocess.run(
                    [cxx] + compile_flags + ["-o", str(binary), str(tmpdir / "correct.cpp")],
                    check=True, capture_output=True, timeout=30,
                )
                break
            except (subprocess.CalledProcessError, FileNotFoundError):
                continue
        else:
            print(f"  Compile error: no working compiler found", file=sys.stderr)
            return 0

        # 各 example に対して correct を実行して .out を生成
        count = 0
        for i, in_name in enumerate(sorted(in_files)):
            in_path = tmpdir / in_name
            out_name = f"sample_{i:02d}"

            try:
                result = subprocess.run(
                    [str(binary)],
                    stdin=open(in_path),
                    capture_output=True, timeout=30,
                )
                if result.returncode != 0:
                    print(f"  {in_name}: correct.cpp returned {result.returncode}", file=sys.stderr)
                    continue

                in_data = in_path.read_text()
                out_data = result.stdout.decode()

                (out_dir / f"{out_name}.in").write_text(in_data)
                (out_dir / f"{out_name}.out").write_text(out_data)
                count += 1
                print(f"  {out_name}: {len(in_data.splitlines())} lines in, {len(out_data.splitlines())} lines out")
            except subprocess.TimeoutExpired:
                print(f"  {in_name}: timeout", file=sys.stderr)
            except Exception as e:
                print(f"  {in_name}: {e}", file=sys.stderr)

        return count


# ============================================================
# AOJ
# ============================================================

def fetch_aoj(url: str, out_dir: Path) -> int:
    """AOJ のサンプルをサンプル専用 API から取得"""
    m = re.search(r"/problems/(\w+)", url) or re.search(r"/(\w+)$", url)
    if not m:
        print(f"  Cannot parse problem ID from {url}", file=sys.stderr)
        return 0
    problem_id = m.group(1)

    api_url = f"https://judgedat.u-aizu.ac.jp/testcases/samples/{problem_id}"
    try:
        with urllib.request.urlopen(api_url, timeout=30) as resp:
            samples = json.loads(resp.read())
    except Exception as e:
        print(f"  AOJ samples API error: {e}", file=sys.stderr)
        return 0

    if not samples:
        print(f"  No samples found for {problem_id}", file=sys.stderr)
        return 0

    count = 0
    for i, tc in enumerate(samples):
        in_data = tc.get("in", "")
        out_data = tc.get("out", "")
        (out_dir / f"sample_{i:02d}.in").write_text(in_data)
        (out_dir / f"sample_{i:02d}.out").write_text(out_data)
        count += 1
        print(f"  sample_{i:02d}: {len(in_data.splitlines())} lines in, {len(out_data.splitlines())} lines out")

    return count


# ============================================================
# yukicoder
# ============================================================

def fetch_yukicoder(url: str, out_dir: Path) -> int:
    """yukicoder のサンプルを問題ページの HTML から取得"""
    m = re.search(r"/problems/no/(\d+)", url)
    if not m:
        print(f"  Cannot parse problem number from {url}", file=sys.stderr)
        return 0
    problem_no = m.group(1)

    page_url = f"https://yukicoder.me/problems/no/{problem_no}"
    req = urllib.request.Request(page_url)
    req.add_header("User-Agent", "Mozilla/5.0")

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            page_html = resp.read().decode()
    except Exception as e:
        print(f"  yukicoder fetch error: {e}", file=sys.stderr)
        return 0

    # id="sample_in_N" / id="sample_out_N" パターン
    inputs = re.findall(r'id=["\']?sample_in_(\d+)["\']?[^>]*>.*?<pre[^>]*>(.*?)</pre>', page_html, re.DOTALL)
    outputs = re.findall(r'id=["\']?sample_out_(\d+)["\']?[^>]*>.*?<pre[^>]*>(.*?)</pre>', page_html, re.DOTALL)

    if not inputs:
        inputs = re.findall(r'id=["\']?sample_in_(\d+)["\']?[^>]*><pre[^>]*>(.*?)</pre>', page_html, re.DOTALL)
        outputs = re.findall(r'id=["\']?sample_out_(\d+)["\']?[^>]*><pre[^>]*>(.*?)</pre>', page_html, re.DOTALL)

    if not inputs:
        print(f"  Could not parse samples from yukicoder page", file=sys.stderr)
        return 0

    out_map = {num: text for num, text in outputs}

    count = 0
    for num, in_html in inputs:
        if num not in out_map:
            continue
        in_text = html_mod.unescape(re.sub(r'<[^>]+>', '', in_html)).strip() + "\n"
        out_text = html_mod.unescape(re.sub(r'<[^>]+>', '', out_map[num])).strip() + "\n"
        (out_dir / f"sample_{count:02d}.in").write_text(in_text)
        (out_dir / f"sample_{count:02d}.out").write_text(out_text)
        count += 1
        print(f"  sample_{count - 1:02d}: {len(in_text.splitlines())} lines in, {len(out_text.splitlines())} lines out")

    return count


# ============================================================
# メイン
# ============================================================

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <problem-dir-or-cpp-file>", file=sys.stderr)
        sys.exit(1)

    target = Path(sys.argv[1])

    # cpp ファイルが指定された場合はディレクトリに変換
    if target.is_file() and target.suffix == ".cpp":
        target = target.parent

    if not target.is_dir():
        print(f"Error: {target} is not a directory", file=sys.stderr)
        sys.exit(1)

    # problem.toml から URL を取得
    config = parse_problem_toml(target)
    url = config.get("url")
    if not url:
        print(f"Error: No url found in {target}/problem.toml", file=sys.stderr)
        sys.exit(1)

    print(f"Problem URL: {url}")

    # testcases/ に保存
    out_dir = target / "testcases"
    out_dir.mkdir(parents=True, exist_ok=True)

    if "judge.yosupo.jp" in url:
        count = fetch_yosupo(url, out_dir)
    elif "onlinejudge.u-aizu.ac.jp" in url or "aoj" in url.lower():
        count = fetch_aoj(url, out_dir)
    elif "yukicoder.me" in url:
        count = fetch_yukicoder(url, out_dir)
    else:
        print(f"Unsupported judge: {url}", file=sys.stderr)
        sys.exit(1)

    if count == 0:
        print("No samples fetched!", file=sys.stderr)
        sys.exit(1)

    print(f"\nFetched {count} sample(s) -> {out_dir}/")


if __name__ == "__main__":
    main()
