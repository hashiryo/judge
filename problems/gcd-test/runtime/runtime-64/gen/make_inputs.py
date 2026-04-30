#!/usr/bin/env python3

from __future__ import annotations

import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTCASES_DIR = ROOT / "testcases"

U64_MAX = (1 << 64) - 1


def write_case(name: str, n: int, abits: int, sa: int, sb: int) -> None:
    TESTCASES_DIR.mkdir(parents=True, exist_ok=True)
    assert 1 <= abits <= 64, f"abits out of range: {abits}"
    path = TESTCASES_DIR / f"{name}.in"
    path.write_text(f"{n} {abits} {sa & U64_MAX} {sb & U64_MAX}\n")


def add_handmade_cases() -> None:
    write_case("small_00", 0, 64, 1, 1)
    write_case("small_01", 1, 64, 1, 1)
    write_case("small_02", 5, 64, 1, 1)
    # 値域の境界: 小さい / 中央 / 最大ビット
    write_case("edge_small", 100_000, 16, 1, 1)
    write_case("edge_mid", 100_000, 32, 1, 1)
    write_case("edge_large", 100_000, 64, 1, 1)


def add_random_cases() -> None:
    rng = random.Random(32)
    configs = [
        ("rand_00", 100_000, 64),
        ("rand_01", 1_000_000, 64),
        ("mid_00", 5_000_000, 64),
        ("mid_01", 5_000_000, 32),
        ("heavy_00", 10_000_000, 64),
        ("heavy_01", 10_000_000, 64),
    ]
    for name, n, abits in configs:
        sa = rng.randrange(U64_MAX + 1)
        sb = rng.randrange(U64_MAX + 1)
        write_case(name, n, abits, sa, sb)


def main() -> None:
    for path in TESTCASES_DIR.glob("*.in"):
        path.unlink()
    for path in TESTCASES_DIR.glob("*.out"):
        path.unlink()

    add_handmade_cases()
    add_random_cases()


if __name__ == "__main__":
    main()
