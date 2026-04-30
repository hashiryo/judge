#!/usr/bin/env python3

from __future__ import annotations

import random
from pathlib import Path

MIN_MOD = (1 << 40) + 1
MAX_MOD = (1 << 62) - 1
ROOT = Path(__file__).resolve().parent.parent
TESTCASES_DIR = ROOT / "testcases"

U64_MAX = (1 << 64) - 1


def normalize_mod(mod: int) -> int:
    mod = max(mod, MIN_MOD)
    mod = min(mod, MAX_MOD)
    if mod % 2 == 0:
        mod += 1 if mod < MAX_MOD else -1
    return mod


def write_case(name: str, n: int, mod: int, bbits: int, sa: int, sb: int) -> None:
    TESTCASES_DIR.mkdir(parents=True, exist_ok=True)
    mod = normalize_mod(mod)
    assert 1 <= bbits <= 64, f"bbits out of range: {bbits}"
    path = TESTCASES_DIR / f"{name}.in"
    path.write_text(f"{n} {mod} {bbits} {sa & U64_MAX} {sb & U64_MAX}\n")


def add_handmade_cases() -> None:
    write_case("small_00", 0, MIN_MOD, 60, 1, 1)
    write_case("small_01", 1, MIN_MOD, 60, 1, 1)
    write_case("small_02", 5, MIN_MOD, 60, 1, 1)
    # bbits の境界
    write_case("edge_b_min", 100_000, MAX_MOD, 1, 1, 1)
    write_case("edge_b_short", 100_000, MAX_MOD, 32, 1, 1)
    write_case("edge_b_max", 100_000, MAX_MOD, 64, 1, 1)
    # mod の境界
    write_case("edge_mod_min", 100_000, MIN_MOD, 60, 1, 1)
    write_case("edge_mod_max", 100_000, MAX_MOD, 60, 1, 1)


def sample_mod(rng: random.Random) -> int:
    mod = rng.randrange(MIN_MOD, MAX_MOD + 1)
    return normalize_mod(mod)


def add_random_cases() -> None:
    rng = random.Random(32)
    # b は 60 ビット固定で、内部 mul は約 60n 回。
    configs = [
        ("rand_00", 10_000, 60),
        ("rand_01", 100_000, 60),
        ("mid_00", 1_000_000, 60),
        ("mid_01", 2_000_000, 60),
        ("heavy_00", 5_000_000, 60),
        ("heavy_01", 5_000_000, 60),
    ]
    for name, n, bbits in configs:
        mod = sample_mod(rng)
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
