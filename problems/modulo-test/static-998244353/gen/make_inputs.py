#!/usr/bin/env python3

from __future__ import annotations

import random
from pathlib import Path

MOD = 998244353
ROOT = Path(__file__).resolve().parent.parent
TESTCASES_DIR = ROOT / "testcases"


def write_case(name: str, n: int, s: int, a: int, b: int, c: int, d: int) -> None:
    TESTCASES_DIR.mkdir(parents=True, exist_ok=True)
    path = TESTCASES_DIR / f"{name}.in"
    path.write_text(f"{n} {s} {a} {b} {c} {d}\n")


def add_handmade_cases() -> None:
    write_case("small_00", 0, 0, 0, 0, 0, 0)
    write_case("small_01", 1, 1, 0, 0, 0, 0)
    write_case("small_02", 5, 1, 1, 1, 1, 1)
    write_case("edge_00", 20, MOD - 1, MOD - 1, 0, MOD - 1, 0)
    write_case("edge_01", 20, MOD - 1, 1, MOD - 1, 1, MOD - 1)
    write_case("edge_02", 100, 123456789, 0, MOD - 1, 0, MOD - 1)


def add_random_cases() -> None:
    rng = random.Random(998244353)
    configs = [
        ("rand_00", 10_000),
        ("rand_01", 100_000),
        ("mid_00", 10_000_000),
        ("mid_01", 30_000_000),
        ("heavy_00", 100_000_000),
        ("heavy_01", 200_000_000),
    ]
    for name, n in configs:
        s = rng.randrange(MOD)
        a = rng.randrange(MOD)
        b = rng.randrange(MOD)
        c = rng.randrange(MOD)
        d = rng.randrange(MOD)
        write_case(name, n, s, a, b, c, d)


def main() -> None:
    for path in TESTCASES_DIR.glob("*.in"):
        path.unlink()
    for path in TESTCASES_DIR.glob("*.out"):
        path.unlink()

    add_handmade_cases()
    add_random_cases()


if __name__ == "__main__":
    main()
