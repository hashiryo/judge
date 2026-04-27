#!/usr/bin/env python3

from __future__ import annotations

import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTCASES_DIR = ROOT / "testcases"

U64_MAX = (1 << 64) - 1
B_MIN = 1
B_MAX = (1 << 63) - 1


def write_case(name: str, n: int, b: int, sa: int, ma: int, mb: int) -> None:
    TESTCASES_DIR.mkdir(parents=True, exist_ok=True)
    assert B_MIN <= b <= B_MAX, f"b out of range: {b}"
    path = TESTCASES_DIR / f"{name}.in"
    path.write_text(f"{n} {b} {sa & U64_MAX} {ma & U64_MAX} {mb & U64_MAX}\n")


def add_handmade_cases() -> None:
    write_case("small_00", 0, B_MIN, 0, 0, 0)
    write_case("small_01", 1, B_MIN, 1, 1, 0)
    write_case("small_02", 5, B_MIN, 1, 6364136223846793005, 1442695040888963407)
    write_case("edge_b_min", 1_000_000, B_MIN, 1, 6364136223846793005, 1442695040888963407)
    write_case("edge_b_max", 1_000_000, B_MAX, 1, 6364136223846793005, 1442695040888963407)
    write_case("edge_b_pow2", 1_000_000, 1 << 40, 1, 6364136223846793005, 1442695040888963407)


def add_random_cases() -> None:
    rng = random.Random(32)
    configs = [
        ("rand_00", 10_000),
        ("rand_01", 100_000),
        ("mid_00", 10_000_000),
        ("mid_01", 30_000_000),
        ("heavy_00", 100_000_000),
        ("heavy_01", 200_000_000),
        ("heavy_02", 200_000_000),
        ("heavy_03", 200_000_000),
    ]
    for name, n in configs:
        b = rng.randrange(B_MIN, B_MAX + 1)
        sa = rng.randrange(U64_MAX + 1)
        ma = rng.randrange(U64_MAX + 1)
        mb = rng.randrange(U64_MAX + 1)
        write_case(name, n, b, sa, ma, mb)


def main() -> None:
    for path in TESTCASES_DIR.glob("*.in"):
        path.unlink()
    for path in TESTCASES_DIR.glob("*.out"):
        path.unlink()

    add_handmade_cases()
    add_random_cases()


if __name__ == "__main__":
    main()
