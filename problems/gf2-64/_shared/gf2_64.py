"""GF(2^64) ≅ GF(2)[x] / (x^64 + x^4 + x^3 + x + 1) の基本演算 (pure Python)。

各 problem の testcases/gen.py から `import gf2_64` する想定。
golden として遅くてよいので素朴実装。
"""
from __future__ import annotations

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
