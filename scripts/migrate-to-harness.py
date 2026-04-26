#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""
旧形式 sol_*.cpp を新形式 base.cpp + algos/*.hpp に変換するワンショット移行スクリプト。

対象は modulo-test 配下のみ (UnionFind 系のような独自構造のものは手動移行)。

使い方:
    python3 scripts/migrate-to-harness.py problems/modulo-test/static/static-1e9+7
    python3 scripts/migrate-to-harness.py --all
"""
import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

COMMON_HPP = """\
#pragma once
// algos 共通: typedef とよく使うヘッダ。
// ファイル名が _ 始まりなので detect-changed.py の提出対象列挙からは除外される。
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;
"""


def make_base_cpp(int_type: str, mod_expr: str | None) -> str:
    """static (mod_expr 指定) と runtime (mod_expr=None) の base.cpp を生成。"""
    if mod_expr is not None:
        mod_decl = f"    constexpr {int_type} MOD = {mod_expr};"
        cin_line = f"    cin >> n_ >> state_ >> a_ >> b_ >> c_ >> d_;"
        var_decl = f"    {int_type} n_, state_, a_, b_, c_, d_;"
        mp_decl = "    const MP mp(MOD);"
    else:
        mod_decl = ""
        cin_line = f"    cin >> n_ >> mod_ >> state_ >> a_ >> b_ >> c_ >> d_;"
        var_decl = f"    {int_type} n_, mod_, state_, a_, b_, c_, d_;"
        mp_decl = "    const MP mp(mod_);"

    pre_loop = "\n".join(filter(None, [mod_decl, var_decl, cin_line]))

    return f"""\
// harness: 各 algos/*.hpp が定義する struct MP を使って計測する
#include "algos/_common.hpp"

// CI では -DALGO_HPP="\\"algos/xxx.hpp\\"" で上書きされる。
// ここでのデフォルトは IDE で base.cpp 単独表示時に補完を効かせるため。
#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {{
    cin.tie(0);
    ios::sync_with_stdio(false);
{pre_loop}

    // REPEAT は wall time を REPEAT 倍するので TLE と相談して決める。
    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    u64 result_out = 0;

    for (int rep = 0; rep < REPEAT; ++rep) {{
{mp_decl}
        auto state = mp.set(state_);
        auto a = mp.set(a_);
        auto b = mp.set(b_);
        auto c = mp.set(c_);
        auto d = mp.set(d_);
        auto t0 = chrono::steady_clock::now();
        for ({int_type} i = 0; i < n_; ++i) {{
            state = mp.mul(mp.plus(mp.mul(state, a), b), mp.plus(mp.mul(state, c), d));
        }}
        auto t1 = chrono::steady_clock::now();
        result_out = mp.get(state);
        auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        if (ns < best_ns) best_ns = ns;
    }}

    fprintf(stderr, "ALGO_TIME_NS=%llu\\n", (unsigned long long)best_ns);
    cout << result_out << '\\n';
    return 0;
}}
"""


# 旧形式ヘッダで除去する行 (typedef / boilerplate)
HEADER_DROP_RE = re.compile(
    r"^\s*(?:"
    r"#include\s*<bits/stdc\+\+\.h>"
    r"|using\s+namespace\s+std\s*;"
    r"|using\s+(?:u8|u32|u64|u128|i64)\s*=\s*[^;]+;"
    r")\s*$"
)


def convert_sol(content: str) -> tuple[str, dict]:
    """
    sol_*.cpp → algos/*.hpp の本文を返す。

    戻り値: (hpp_body, info) — info = {"has_mod_input": bool, "int_type": str}
    """
    lines = content.splitlines()

    # main の開始行を探す
    main_idx = None
    for i, line in enumerate(lines):
        if re.match(r"^\s*(?:signed|int)\s+main\s*\(", line):
            main_idx = i
            break
    if main_idx is None:
        raise ValueError("main() が見つからない")

    pre = lines[:main_idx]
    body = "\n".join(lines[main_idx:])

    # main から mp 宣言の型を抽出
    mp_match = re.search(r"(?:constexpr\s+)?(\w+)\s+mp\s*\(", body)
    if not mp_match:
        raise ValueError("mp(...) 宣言が見つからない")
    mp_type = mp_match.group(1)

    # 入力パターン (runtime かどうか)
    has_mod_input = bool(re.search(r"cin\s*>>\s*n\s*>>\s*mod\b", body))

    # 整数型 (u32 / u64) — main 内の `u32 n,` / `u64 n,` から
    type_match = re.search(r"\b(u32|u64)\s+n\s*,", body)
    int_type = type_match.group(1) if type_match else "u32"

    # ヘッダ行を落とす
    kept = []
    leading_blank = True
    for line in pre:
        if HEADER_DROP_RE.match(line):
            continue
        # 先頭の空行は落とす
        if leading_blank and not line.strip():
            continue
        leading_blank = False
        kept.append(line)
    # 末尾の空行も落とす
    while kept and not kept[-1].strip():
        kept.pop()

    pre_text = "\n".join(kept)

    # mp_type → MP に置換 (単語境界で)
    pre_text = re.sub(rf"\b{re.escape(mp_type)}\b", "MP", pre_text)

    hpp = f"""#pragma once
#include "_common.hpp"
{pre_text}
"""
    return hpp, {"has_mod_input": has_mod_input, "int_type": int_type}


def algo_name_from_filename(name: str) -> str:
    """sol_xxx.cpp → xxx (ファイル名として安全な形)"""
    stem = Path(name).stem
    if stem.startswith("sol_"):
        stem = stem[4:]
    # 等号などをアンダースコアに置換 (ファイル名としては OK だが念のため)
    return stem


# 問題ディレクトリ → 静的 MOD 式 (None なら runtime)
PROBLEM_MOD_TABLE: dict[str, str | None] = {
    "problems/modulo-test/static/static-998244353": "998244353",
    "problems/modulo-test/static/static-1e9+7": "1000000007",
    "problems/modulo-test/static/static-2^31-1": "(1u << 31) - 1",
    "problems/modulo-test/runtime/runtime-30": None,
    "problems/modulo-test/runtime/runtime-31": None,
    "problems/modulo-test/runtime/runtime-32": None,
    "problems/modulo-test/runtime/runtime-40": None,
}


def migrate_problem(prob_dir: Path, dry_run: bool = False) -> None:
    rel = str(prob_dir.relative_to(ROOT))
    if rel not in PROBLEM_MOD_TABLE:
        print(f"  skip (未登録): {rel}", file=sys.stderr)
        return

    mod_expr = PROBLEM_MOD_TABLE[rel]
    sol_files = sorted(prob_dir.glob("sol_*.cpp"))
    if not sol_files:
        print(f"  skip (sol_*.cpp なし): {rel}", file=sys.stderr)
        return

    algos_dir = prob_dir / "algos"
    int_types: set[str] = set()
    has_mod_input_set: set[bool] = set()
    converted: list[tuple[Path, str]] = []

    for sol in sol_files:
        try:
            hpp_body, info = convert_sol(sol.read_text())
        except Exception as e:
            print(f"  ERROR ({sol.name}): {e}", file=sys.stderr)
            continue
        algo = algo_name_from_filename(sol.name)
        algo_path = algos_dir / f"{algo}.hpp"
        converted.append((algo_path, hpp_body))
        int_types.add(info["int_type"])
        has_mod_input_set.add(info["has_mod_input"])

    if not converted:
        return

    # base.cpp の型: 一番広い型に合わせる
    base_int_type = "u64" if "u64" in int_types else "u32"
    # 入力の mod 有無は問題ディレクトリで一意のはず
    runtime_inputs = (mod_expr is None)
    if has_mod_input_set != {runtime_inputs}:
        print(
            f"  WARN: {rel} の sol_*.cpp が runtime/static で混在しているように見える "
            f"(has_mod_input={has_mod_input_set}, expected={runtime_inputs})",
            file=sys.stderr,
        )

    base_cpp = make_base_cpp(base_int_type, mod_expr)
    common_hpp = COMMON_HPP

    # 出力
    print(f"== {rel} ({len(converted)} algos, type={base_int_type}, mod={'runtime' if runtime_inputs else mod_expr}) ==")
    if dry_run:
        for path, _ in converted:
            print(f"  would write {path.relative_to(ROOT)}")
        print(f"  would write {(prob_dir / 'base.cpp').relative_to(ROOT)}")
        print(f"  would write {(algos_dir / '_common.hpp').relative_to(ROOT)}")
        return

    algos_dir.mkdir(parents=True, exist_ok=True)
    common_path = algos_dir / "_common.hpp"
    if common_path.exists():
        print(f"  keep {common_path.relative_to(ROOT)} (既存)")
    else:
        common_path.write_text(common_hpp)
        print(f"  wrote {common_path.relative_to(ROOT)}")
    for path, body in converted:
        if path.exists():
            print(f"  keep {path.relative_to(ROOT)} (既存)")
        else:
            path.write_text(body)
            print(f"  wrote {path.relative_to(ROOT)}")
    base_path = prob_dir / "base.cpp"
    if base_path.exists():
        print(f"  keep {base_path.relative_to(ROOT)} (既存)")
    else:
        base_path.write_text(base_cpp)
        print(f"  wrote {base_path.relative_to(ROOT)}")

    # 旧 sol_*.cpp を削除
    for sol in sol_files:
        sol.unlink()
        print(f"  removed {sol.relative_to(ROOT)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs="*", help="問題ディレクトリ (relative to repo root)")
    parser.add_argument("--all", action="store_true", help="登録済み全問題を移行")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.all:
        targets = [ROOT / p for p in PROBLEM_MOD_TABLE]
    else:
        if not args.targets:
            parser.error("ディレクトリか --all を指定してください")
        targets = [Path(t).resolve() for t in args.targets]

    for t in targets:
        if not t.is_dir():
            print(f"  not a dir: {t}", file=sys.stderr)
            continue
        migrate_problem(t, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
