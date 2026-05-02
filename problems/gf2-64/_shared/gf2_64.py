"""GF(2^64) ≅ GF(2)[x] / (x^64 + x^4 + x^3 + x + 1) の基本演算。

各 problem の testcases/gen.py から `import gf2_64` する想定。
golden 単体は素朴実装 (gf_mul / gf_pow / gf_inv / gf_sqrt) だが、
T 個一括の batch 呼び出しは scalar C++ helper (gen_helper.cpp) を使う高速版。
T=1M でも数秒で済むので CI gen 時間が現実的になる。
"""
from __future__ import annotations
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# 既約多項式 P(x) = x^64 + x^4 + x^3 + x + 1。下位 64 bit + 暗黙の x^64。
P_LOW = (1 << 4) | (1 << 3) | (1 << 1) | 1
MASK64 = (1 << 64) - 1

# 2^64 - 1 = 18446744073709551615 の素因数 (Mersenne 風 Fermat 数の積)。
# 原始元判定 (g^((2^64-1)/p) != 1 for all p) に使う。
ORDER_PRIMES = (3, 5, 17, 257, 641, 65537, 6700417)
ORDER = MASK64  # 乗法群の位数 = 2^64 - 1


def gf_clmul(a: int, b: int) -> tuple[int, int]:
    """carryless multiply: 64 × 64 → 128 (lo, hi)。"""
    lo = 0
    hi = 0
    for i in range(64):
        if (b >> i) & 1:
            lo ^= (a << i) & MASK64
            if i:
                hi ^= a >> (64 - i)
    return lo, hi


def gf_reduce(lo: int, hi: int) -> int:
    """(lo, hi) を P(x) で reduce。素朴長除法。"""
    # bit i (i ≥ 64) を順に消す: hi の bit i → lo に P_LOW * x^i を xor
    for i in range(63, -1, -1):
        if (hi >> i) & 1:
            hi ^= 1 << i
            shifted_lo = (P_LOW << i) & MASK64
            shifted_hi = P_LOW >> (64 - i) if i > 0 else 0
            lo ^= shifted_lo
            hi ^= shifted_hi
    return lo


def gf_mul(a: int, b: int) -> int:
    return gf_reduce(*gf_clmul(a, b))


def gf_pow(a: int, e: int) -> int:
    res = 1
    while e:
        if e & 1:
            res = gf_mul(res, a)
        a = gf_mul(a, a)
        e >>= 1
    return res


def gf_inv(a: int) -> int:
    """Fermat: a^(2^64 - 2)。a != 0 前提。"""
    return gf_pow(a, MASK64 - 1)


def gf_sqrt(a: int) -> int:
    """sqrt(a) = a^(2^63)。"""
    return gf_pow(a, 1 << 63)


def is_primitive(g: int) -> bool:
    """g が乗法群の原始元か (g^((2^64-1)/p) != 1 for each prime divisor p)。"""
    if g == 0:
        return False
    for p in ORDER_PRIMES:
        if gf_pow(g, ORDER // p) == 1:
            return False
    return True


def find_primitive(start: int = 2) -> int:
    """start から順に試して最初に見つかった原始元を返す。"""
    g = start
    while True:
        if is_primitive(g):
            return g
        g += 1


# log 用に 1 度だけ確定する固定の原始元。問題定義の一部としてハードコード。
# (問題本文・checker・gen 全てで同じ値を使う)
PRIMITIVE_G = 2  # 仮置き。下で is_primitive 確認後、必要なら更新。
if not is_primitive(PRIMITIVE_G):
    PRIMITIVE_G = find_primitive()


# ============================================================
# C++ helper (gen_helper.cpp) を介する batch 演算。
# pure Python の gf_mul は 64-loop × 1M 個 = ~12 秒かかるので、
# 同等の scalar C++ 実装を 1 度コンパイルして batch でまとめて流す。
# ============================================================

_HERE = Path(__file__).resolve().parent
_HELPER_SRC = _HERE / "gen_helper.cpp"
_HELPER_BIN_CACHE: Path | None = None


def _ensure_helper_compiled() -> Path:
    """gen_helper を一度だけコンパイル。/tmp 上に置いて以降は再利用。"""
    global _HELPER_BIN_CACHE
    if _HELPER_BIN_CACHE is not None and _HELPER_BIN_CACHE.exists():
        return _HELPER_BIN_CACHE

    # ソースの mtime を hash に混ぜて、ソース更新時は再コンパイル。
    src_mtime = int(_HELPER_SRC.stat().st_mtime)
    bin_path = Path(tempfile.gettempdir()) / f"gf2_64_gen_helper.{src_mtime}.bin"

    if not bin_path.exists():
        cxx = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
        if cxx is None:
            raise RuntimeError("no C++ compiler found (c++/g++/clang++)")
        print(f"  [gf2_64] compiling helper -> {bin_path.name}", file=sys.stderr)
        subprocess.run(
            [cxx, "-std=c++17", "-O3", "-o", str(bin_path), str(_HELPER_SRC)],
            check=True,
        )

    _HELPER_BIN_CACHE = bin_path
    return bin_path


def _run_helper(op: str, lines: list[str]) -> list[str]:
    """helper を起動し、stdin に op + 入力行を流し、stdout から T 行受け取る。"""
    bin_path = _ensure_helper_compiled()
    payload = f"{op} {len(lines)}\n" + "\n".join(lines) + "\n"
    res = subprocess.run(
        [str(bin_path)], input=payload, capture_output=True, text=True, check=True,
    )
    out = res.stdout.splitlines()
    if len(out) != len(lines):
        raise RuntimeError(f"helper {op}: expected {len(lines)} lines, got {len(out)}")
    return out


def batch_mul(pairs: list[tuple[int, int]]) -> list[int]:
    out = _run_helper("mul", [f"{a} {b}" for a, b in pairs])
    return [int(x) for x in out]


def batch_div(pairs: list[tuple[int, int]]) -> list[int]:
    """(a, b) → a/b。b != 0 前提。"""
    out = _run_helper("div", [f"{a} {b}" for a, b in pairs])
    return [int(x) for x in out]


def batch_pow(pairs: list[tuple[int, int]]) -> list[int]:
    out = _run_helper("pow", [f"{a} {e}" for a, e in pairs])
    return [int(x) for x in out]


def batch_sqrt(values: list[int]) -> list[int]:
    out = _run_helper("sqrt", [str(a) for a in values])
    return [int(x) for x in out]


def batch_inv(values: list[int]) -> list[int]:
    out = _run_helper("inv", [str(a) for a in values])
    return [int(x) for x in out]


def batch_log_inv(ks: list[int]) -> list[int]:
    """k → g^k (g = 2)。log 問題のテストケース生成で「正解 k から x を作る」のに使う。"""
    out = _run_helper("log_inv", [str(k) for k in ks])
    return [int(x) for x in out]
