#!/usr/bin/env python3

from __future__ import annotations

import random
from pathlib import Path

MIN_MOD = (1 << 31) + 1
MAX_MOD = (1 << 32) - 1
ROOT = Path(__file__).resolve().parent.parent
TESTCASES_DIR = ROOT / "testcases"


def normalize_mod(mod: int) -> int:
    mod = max(mod, MIN_MOD)
    mod = min(mod, MAX_MOD)
    if mod % 2 == 0:
        mod += 1 if mod < MAX_MOD else -1
    return mod


def write_case(name: str, n: int, mod: int, s: int, a: int, b: int, c: int, d: int) -> None:
    TESTCASES_DIR.mkdir(parents=True, exist_ok=True)
    mod = normalize_mod(mod)
    path = TESTCASES_DIR / f"{name}.in"
    path.write_text(f"{n} {mod} {s % mod} {a % mod} {b % mod} {c % mod} {d % mod}\n")


def add_handmade_cases() -> None:
    write_case("small_00", 0, MIN_MOD, 0, 0, 0, 0, 0)
    write_case("small_01", 1, MIN_MOD, 1, 0, 0, 0, 0)
    write_case("small_02", 5, MIN_MOD, 1, 1, 1, 1, 1)
    write_case("edge_00", 20, MIN_MOD, MIN_MOD - 1, MIN_MOD - 1, 0, MIN_MOD - 1, 0)
    write_case("edge_01", 20, (1 << 31) - 1, (1 << 31) - 2, 1, (1 << 31) - 2, 1, (1 << 31) - 2)
    write_case("edge_02", 100, MAX_MOD, 123456789, 0, MAX_MOD - 1, 0, MAX_MOD - 1)


def sample_mod(rng: random.Random) -> int:
    mod = rng.randrange(MIN_MOD, MAX_MOD + 1)
    return normalize_mod(mod)


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
        mod = sample_mod(rng)
        s = rng.randrange(mod)
        a = rng.randrange(mod)
        b = rng.randrange(mod)
        c = rng.randrange(mod)
        d = rng.randrange(mod)
        write_case(name, n, mod, s, a, b, c, d)


def main() -> None:
    for path in TESTCASES_DIR.glob("*.in"):
        path.unlink()
    for path in TESTCASES_DIR.glob("*.out"):
        path.unlink()

    add_handmade_cases()
    add_random_cases()


if __name__ == "__main__":
    main()
