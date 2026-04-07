#!/usr/bin/env python3
"""
問題ディレクトリの .cpp ファイルから judge URL を読み取り、テストケースをダウンロードする。
自作テストケース (testcases/) も検出してキャッシュにコピーする。

Library リポジトリの download-testcases.py をベースに judge 用に改変。
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TC_CACHE_DIR = ROOT / ".cache" / "testcases"
CUSTOM_TC_DIR = ROOT / ".cache" / "custom-testcases"

YUKICODER_TOKEN = os.environ.get("YUKICODER_TOKEN", "")


def url_to_md5(url: str) -> str:
    return hashlib.md5(url.encode()).hexdigest()


def url_to_cache_dir(url: str) -> Path:
    return TC_CACHE_DIR / url_to_md5(url)


def parse_judge_url(cpp_file: Path) -> str | None:
    """cpp ファイルの先頭から // judge: URL を読み取る"""
    try:
        with open(cpp_file) as f:
            for line in f:
                line = line.strip()
                if not line.startswith("//"):
                    if line and not line.startswith("#"):
                        break
                    continue
                m = re.match(r"//\s*judge:\s*(\S+)", line)
                if m:
                    return m.group(1)
    except Exception:
        pass
    return None


def collect_problem_urls(problems_json: dict) -> dict[str, list[str]]:
    """問題 JSON から URL → ファイルリストのマッピングを作成"""
    urls: dict[str, list[str]] = {}
    for problem in problems_json.get("problems", []):
        problem_dir = ROOT / problem["dir"]
        for filename in problem["files"]:
            cpp_file = problem_dir / filename
            url = parse_judge_url(cpp_file)
            if url:
                if url not in urls:
                    urls[url] = []
                urls[url].append(str(cpp_file.relative_to(ROOT)))
    return urls


def collect_custom_testcases(problems_json: dict) -> list[tuple[str, Path]]:
    """自作テストケースを持つ問題ディレクトリを検出"""
    result = []
    for problem in problems_json.get("problems", []):
        problem_dir = ROOT / problem["dir"]
        tc_dir = problem_dir / "testcases"
        if tc_dir.is_dir():
            in_files = list(tc_dir.glob("*.in"))
            if in_files:
                result.append((problem["dir"], tc_dir))
    return result


def setup_custom_testcases(problem_dir_name: str, tc_dir: Path):
    """自作テストケースをキャッシュディレクトリにコピー"""
    # 自作テストケースは問題ディレクトリ名をキーにキャッシュ
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

    # checker があればコピー
    checker = tc_dir / "checker.cpp"
    if checker.exists():
        shutil.copy2(checker, cache_dir / "checker.cpp")

    print(f"  Custom testcases for {problem_dir_name}: {count} cases")


def is_cached(url: str) -> bool:
    cache_dir = url_to_cache_dir(url)
    return cache_dir.exists() and any(cache_dir.iterdir())


# ============================================================
# AOJ
# ============================================================

def download_aoj(url: str) -> bool:
    cache_dir = url_to_cache_dir(url)
    m = re.search(r"/problems/(\w+)", url) or re.search(r"/(\w+)$", url)
    if not m:
        return False
    problem_id = m.group(1)

    try:
        header_url = f"https://judgedat.u-aizu.ac.jp/testcases/{problem_id}/header"
        with urllib.request.urlopen(header_url, timeout=30) as resp:
            headers = json.loads(resp.read())
    except Exception as e:
        print(f"    AOJ {problem_id}: {e}", flush=True)
        return False

    tmp_dir = Path(str(cache_dir) + ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    expected = [tc for tc in headers.get("headers", []) if tc.get("serial") is not None]
    count = 0

    for tc in expected:
        serial = tc["serial"]
        name = tc.get("name", f"case{serial}")
        try:
            in_url = f"https://judgedat.u-aizu.ac.jp/testcases/{problem_id}/{serial}/in"
            with urllib.request.urlopen(in_url, timeout=30) as resp:
                in_data = resp.read()
            out_url = f"https://judgedat.u-aizu.ac.jp/testcases/{problem_id}/{serial}/out"
            with urllib.request.urlopen(out_url, timeout=30) as resp:
                out_data = resp.read()
            (tmp_dir / f"{name}.in").write_bytes(in_data)
            (tmp_dir / f"{name}.out").write_bytes(out_data)
            count += 1
        except Exception:
            shutil.rmtree(tmp_dir)
            return False

    if count > 0:
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
        tmp_dir.rename(cache_dir)
        return True
    else:
        shutil.rmtree(tmp_dir)
        return False


# ============================================================
# yosupo judge (Library Checker)
# ============================================================

LIBRARY_CHECKER_DIR = ROOT / ".cache" / "library-checker-problems"


def ensure_library_checker_repo():
    if LIBRARY_CHECKER_DIR.exists():
        print("  Updating library-checker-problems...", flush=True)
        try:
            subprocess.run(
                ["git", "fetch", "--depth=1", "origin", "master"],
                check=True, capture_output=True, timeout=60,
                cwd=LIBRARY_CHECKER_DIR,
            )
            subprocess.run(
                ["git", "reset", "--hard", "origin/master"],
                check=True, capture_output=True, timeout=30,
                cwd=LIBRARY_CHECKER_DIR,
            )
            return True
        except Exception as e:
            print(f"  Update failed: {e}", flush=True)
            return True
    else:
        print("  Cloning library-checker-problems...", flush=True)
        try:
            subprocess.run(
                ["git", "clone", "--depth=1",
                 "https://github.com/yosupo06/library-checker-problems.git",
                 str(LIBRARY_CHECKER_DIR)],
                check=True, capture_output=True, timeout=120,
            )
        except Exception as e:
            print(f"  Clone failed: {e}", flush=True)
            return False

    subprocess.run(
        [sys.executable, "-m", "pip", "install", "-q", "pyyaml", "toml"],
        capture_output=True,
    )
    return True


def download_yosupo(url: str) -> bool:
    cache_dir = url_to_cache_dir(url)
    m = re.search(r"/problem/(\w+)", url)
    if not m:
        return False
    problem_id = m.group(1)

    if not ensure_library_checker_repo():
        return False

    problem_dir = None
    for p in LIBRARY_CHECKER_DIR.rglob("info.toml"):
        if p.parent.name == problem_id:
            problem_dir = p.parent
            break
    if not problem_dir:
        return False

    try:
        subprocess.run(
            ["bash", "-c",
             f"ulimit -s unlimited && {sys.executable} {LIBRARY_CHECKER_DIR / 'generate.py'} {problem_dir / 'info.toml'}"],
            check=True, timeout=300,
            cwd=LIBRARY_CHECKER_DIR,
        )
    except Exception as e:
        print(f"    yosupo generate failed for {problem_id}: {e}", flush=True)
        return False

    in_dir = problem_dir / "in"
    out_dir = problem_dir / "out"
    if not in_dir.exists():
        return False

    tmp_dir = Path(str(cache_dir) + ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for in_file in sorted(in_dir.glob("*.in")):
        out_file = out_dir / in_file.name.replace(".in", ".out")
        if out_file.exists():
            shutil.copy2(in_file, tmp_dir / in_file.name)
            shutil.copy2(out_file, tmp_dir / out_file.name)
            count += 1

    checker_cpp = problem_dir / "checker.cpp"
    if checker_cpp.exists():
        shutil.copy2(checker_cpp, tmp_dir / "checker.cpp")
        testlib_h = LIBRARY_CHECKER_DIR / "common" / "testlib.h"
        if testlib_h.exists():
            shutil.copy2(testlib_h, tmp_dir / "testlib.h")

    if count > 0:
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
        tmp_dir.rename(cache_dir)
        return True
    else:
        shutil.rmtree(tmp_dir)
        return False


# ============================================================
# yukicoder
# ============================================================

def download_yukicoder(url: str) -> bool:
    import zipfile
    import io

    cache_dir = url_to_cache_dir(url)
    if not YUKICODER_TOKEN:
        print(f"    yukicoder: no token set", flush=True)
        return False

    m = re.search(r"/problems/no/(\d+)", url)
    if not m:
        print(f"    yukicoder: cannot parse problem number from {url}", flush=True)
        return False
    problem_no = m.group(1)

    zip_url = f"https://yukicoder.me/problems/no/{problem_no}/testcase.zip"
    req = urllib.request.Request(zip_url)
    req.add_header("Authorization", f"Bearer {YUKICODER_TOKEN}")

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            zip_data = resp.read()
    except Exception as e:
        print(f"    yukicoder #{problem_no}: {e}", flush=True)
        return False

    tmp_dir = Path(str(cache_dir) + ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    try:
        with zipfile.ZipFile(io.BytesIO(zip_data)) as zf:
            inputs = {}
            outputs = {}
            for name in zf.namelist():
                if name.startswith("test_in/") and not name.endswith("/"):
                    case_name = Path(name).stem
                    inputs[case_name] = zf.read(name)
                elif name.startswith("test_out/") and not name.endswith("/"):
                    case_name = Path(name).stem
                    outputs[case_name] = zf.read(name)

            count = 0
            for case_name in inputs:
                if case_name in outputs:
                    (tmp_dir / f"{case_name}.in").write_bytes(inputs[case_name])
                    (tmp_dir / f"{case_name}.out").write_bytes(outputs[case_name])
                    count += 1
    except Exception as e:
        print(f"    yukicoder #{problem_no}: zip extract failed: {e}", flush=True)
        shutil.rmtree(tmp_dir)
        return False

    if count > 0:
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
        tmp_dir.rename(cache_dir)
        return True
    else:
        shutil.rmtree(tmp_dir)
        return False


# ============================================================
# HackerRank
# ============================================================

def download_hackerrank(url: str) -> bool:
    import zipfile
    import io

    cache_dir = url_to_cache_dir(url)
    m = re.search(r"hackerrank\.com/contests/([^/]+)/challenges/([^/]+)", url)
    if m:
        contest, challenge = m.group(1), m.group(2)
    else:
        m = re.search(r"hackerrank\.com/challenges/([^/]+)", url)
        if not m:
            return False
        contest, challenge = "master", m.group(1)

    api_url = f"https://www.hackerrank.com/rest/contests/{contest}/challenges/{challenge}/download_testcases"
    req = urllib.request.Request(api_url)
    req.add_header("User-Agent", "Mozilla/5.0")

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            zip_data = resp.read()
    except Exception as e:
        print(f"    HackerRank {challenge}: {e}", flush=True)
        return False

    tmp_dir = Path(str(cache_dir) + ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    try:
        with zipfile.ZipFile(io.BytesIO(zip_data)) as zf:
            count = 0
            for name in zf.namelist():
                if name.startswith("input/") and name.endswith(".txt"):
                    case_name = Path(name).stem
                    in_data = zf.read(name)
                    out_name = name.replace("input/", "output/").replace("input", "output")
                    if out_name in zf.namelist():
                        out_data = zf.read(out_name)
                        (tmp_dir / f"{case_name}.in").write_bytes(in_data)
                        (tmp_dir / f"{case_name}.out").write_bytes(out_data)
                        count += 1
    except Exception:
        shutil.rmtree(tmp_dir)
        return False

    if count > 0:
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
        tmp_dir.rename(cache_dir)
        return True
    else:
        shutil.rmtree(tmp_dir)
        return False


# ============================================================
# ディスパッチ
# ============================================================

def download_one(url: str) -> bool:
    if "onlinejudge.u-aizu.ac.jp" in url:
        return download_aoj(url)
    elif "judge.yosupo.jp" in url:
        return download_yosupo(url)
    elif "yukicoder.me" in url:
        return download_yukicoder(url)
    elif "hackerrank.com" in url:
        return download_hackerrank(url)
    else:
        return False


def download_one_safe(url: str) -> tuple[str, bool]:
    try:
        return (url, download_one(url))
    except Exception as e:
        print(f"    Error: {url}: {e}", file=sys.stderr)
        return (url, False)


# ============================================================
# メイン
# ============================================================

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--problems-json", required=True, help="detect-changed.py の出力 JSON")
    args = parser.parse_args()

    problems_json = json.loads(args.problems_json)
    urls = collect_problem_urls(problems_json)
    print(f"Found {len(urls)} unique problem URLs")

    # 自作テストケースの処理
    custom_tcs = collect_custom_testcases(problems_json)
    for problem_dir_name, tc_dir in custom_tcs:
        setup_custom_testcases(problem_dir_name, tc_dir)

    if not urls:
        if custom_tcs:
            print("Only custom testcases, no downloads needed")
        else:
            print("No problem URLs found")
        return

    cached = 0
    to_download: list[str] = []

    for url in urls:
        if is_cached(url):
            cached += 1
        else:
            to_download.append(url)

    print(f"  Cached: {cached}")
    print(f"  Need download: {len(to_download)}")

    if not to_download:
        print("All testcases available!")
        return

    yosupo_urls = [u for u in to_download if "judge.yosupo.jp" in u]
    other_urls = [u for u in to_download if "judge.yosupo.jp" not in u]

    downloaded = 0
    failed = 0

    if other_urls:
        print(f"  Downloading {len(other_urls)} problems (parallel)...")
        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = {executor.submit(download_one_safe, url): url for url in other_urls}
            for future in as_completed(futures):
                url, success = future.result()
                if success:
                    downloaded += 1
                    print(f"  [OK]  {url}")
                else:
                    failed += 1
                    print(f"  [NG]  {url}")

    if yosupo_urls:
        if not ensure_library_checker_repo():
            print("  Failed to prepare library-checker-problems, skipping yosupo")
            failed += len(yosupo_urls)
        else:
            print(f"  Generating {len(yosupo_urls)} yosupo problems (parallel)...")
            with ProcessPoolExecutor(max_workers=1) as executor:
                futures = {executor.submit(download_one_safe, url): url for url in yosupo_urls}
                for future in as_completed(futures):
                    url, success = future.result()
                    m = re.search(r"/problem/(\w+)", url)
                    name = m.group(1) if m else url
                    if success:
                        downloaded += 1
                        print(f"  [GEN] {name} ... OK", flush=True)
                    else:
                        failed += 1
                        print(f"  [GEN] {name} ... FAILED", flush=True)

    print(f"\nTestcase download summary:")
    print(f"  Cached:     {cached}")
    print(f"  Downloaded: {downloaded}")
    print(f"  Failed:     {failed}")
    print(f"  Total URLs: {len(urls)}")
    print(f"  Custom:     {len(custom_tcs)} problem(s)")


if __name__ == "__main__":
    main()
