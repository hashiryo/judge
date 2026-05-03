#!/usr/bin/env python3

from __future__ import annotations

import random
from pathlib import Path

# 30-bit prime range. 998244353 (= 119 * 2^23 + 1) を中心に使う。
PRIMES_30BIT = [
    998244353,  # NTT-friendly prime
    1000000007,
    1000000009,
    924844033,
    985661441,
    754974721,
    167772161,
]
DEFAULT_PRIME = 998244353

ROOT = Path(__file__).resolve().parent.parent
TESTCASES_DIR = ROOT / "testcases"

U64_MAX = (1 << 64) - 1


def write_case(name: str, n: int, mod: int, bbits: int, sa: int, sb: int) -> None:
    TESTCASES_DIR.mkdir(parents=True, exist_ok=True)
    assert 1 <= bbits <= 32, f"bbits out of range: {bbits}"
    path = TESTCASES_DIR / f"{name}.in"
    path.write_text(f"{n} {mod} {bbits} {sa & U64_MAX} {sb & U64_MAX}\n")


def add_handmade_cases() -> None:
    write_case("small_00", 0, DEFAULT_PRIME, 30, 1, 1)
    write_case("small_01", 1, DEFAULT_PRIME, 30, 1, 1)
    write_case("small_02", 5, DEFAULT_PRIME, 30, 1, 1)
    # bbits の境界
    write_case("edge_b_min", 100_000, DEFAULT_PRIME, 1, 1, 1)
    write_case("edge_b_short", 100_000, DEFAULT_PRIME, 16, 1, 1)
    write_case("edge_b_med", 100_000, DEFAULT_PRIME, 30, 1, 1)
    write_case("edge_b_max", 100_000, DEFAULT_PRIME, 32, 1, 1)


def add_random_cases() -> None:
    rng = random.Random(32)
    # b は 30 ビット (≈ p-1 と同程度) と 32 ビット (フェルマー削減が効く) の両方。
    configs = [
        ("rand_00", 10_000, 30),
        ("rand_01", 100_000, 30),
        ("mid_00", 1_000_000, 30),
        ("mid_01", 1_000_000, 32),  # b > p-1 の領域 (フェルマー削減で短縮可)
        ("heavy_00", 5_000_000, 30),
        ("heavy_01", 5_000_000, 32),
    ]
    for name, n, bbits in configs:
        mod = rng.choice(PRIMES_30BIT)
        sa = rng.randrange(U64_MAX + 1)
        sb = rng.randrange(U64_MAX + 1)
        write_case(name, n, mod, bbits, sa, sb)


def main() -> None:
    for path in TESTCASES_DIR.glob("*.in"):
        path.unlink()
    for path in TESTCASES_DIR.glob("*.out"):
        path.unlink()

    add_handmade_cases()
    add_random_cases()


if __name__ == "__main__":
    main()
