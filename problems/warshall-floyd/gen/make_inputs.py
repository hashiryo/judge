#!/usr/bin/env python3

from __future__ import annotations

import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTCASES_DIR = ROOT / "testcases"

U64_MAX = (1 << 64) - 1


def write_case(name: str, V: int, W: int, seed: int) -> None:
    TESTCASES_DIR.mkdir(parents=True, exist_ok=True)
    assert 1 <= V, f"V must be >= 1, got {V}"
    assert 1 <= W, f"W must be >= 1, got {W}"
    path = TESTCASES_DIR / f"{name}.in"
    path.write_text(f"{V} {W} {seed & U64_MAX}\n")


def add_handmade_cases() -> None:
    write_case("small_00", 1, 100, 1)
    write_case("small_01", 5, 100, 1)
    write_case("small_02", 10, 100, 1)
    write_case("edge_v8", 8, 100, 1)        # SIMD 1 ベクトル分ぴったり
    write_case("edge_v9", 9, 100, 1)        # パディング 7 列
    write_case("edge_v15", 15, 100, 1)      # パディング 1 列


def add_random_cases() -> None:
    rng = random.Random(32)
    configs = [
        ("rand_00", 50),
        ("rand_01", 100),
        ("mid_00", 256),     # SIMD 整列
        ("mid_01", 300),
        ("heavy_00", 512),   # SIMD 整列、512^3 = 1.34e8
        ("heavy_01", 768),   # 768^3 = 4.5e8
    ]
    for name, V in configs:
        seed = rng.randrange(U64_MAX + 1)
        write_case(name, V, 100, seed)


def main() -> None:
    for path in TESTCASES_DIR.glob("*.in"):
        path.unlink()
    for path in TESTCASES_DIR.glob("*.out"):
        path.unlink()

    add_handmade_cases()
    add_random_cases()


if __name__ == "__main__":
    main()
