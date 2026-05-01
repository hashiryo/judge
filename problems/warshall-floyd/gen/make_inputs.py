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
    # W は i32 dist のオーバフローを避けるため、V * W < WF_INF = 10^9 を満たす範囲。
    # 余裕みて max_path = V * W <= 10^8 程度を狙う。
    configs = [
        # (name, V, W) — W を小さめ / 中間 / 大きめに散らす
        ("rand_w_small", 100, 10),         # 重み小、経路が縮みやすい
        ("rand_w_mid", 100, 10**3),
        ("rand_w_large", 100, 10**5),      # 競プロ AOJ 系想定
        ("mid_w_small", 256, 10),
        ("mid_w_mid", 300, 10**4),
        ("mid_w_large", 256, 10**5),
        ("heavy_00", 512, 10**3),          # 512^3 = 1.34e8 ops
        ("heavy_01", 768, 10**3),          # 768^3 = 4.5e8 ops
        ("heavy_02", 768, 10**5),          # 重みも大きめ
    ]
    for name, V, W in configs:
        seed = rng.randrange(U64_MAX + 1)
        write_case(name, V, W, seed)


def main() -> None:
    for path in TESTCASES_DIR.glob("*.in"):
        path.unlink()
    for path in TESTCASES_DIR.glob("*.out"):
        path.unlink()

    add_handmade_cases()
    add_random_cases()


if __name__ == "__main__":
    main()
