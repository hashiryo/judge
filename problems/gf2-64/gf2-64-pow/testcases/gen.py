#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""GF(2^64) 累乗 a^e のテストケース生成。"""
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "_shared"))
import gf2_64  # noqa: E402

MASK64 = gf2_64.MASK64

CASES = {
    "sample_00":           (10,    "sample"),
    "small_00":            (100,   "small_e"),
    "edge_e_00":           (200,   "edge_e"),
    "random_00":           (10000, "random"),
    "random_large_00":     (100_000, "random"),  # pow は 1 つあたり ~64 mul なので 1e5 で十分
}


def make_pairs(name: str, T: int, kind: str) -> list[tuple[int, int]]:
    rng = random.Random(hash(("pow", name)) & 0xFFFFFFFF)
    if kind == "sample":
        return [
            (0, 0), (0, 1), (1, 5), (2, 0), (2, 1),
            (2, 64), (3, MASK64), (5, 1 << 63),
            (MASK64, MASK64), (0xDEADBEEF, 0x1234567),
        ][:T]
    if kind == "small_e":
        return [(rng.randint(0, MASK64), rng.randint(0, 100)) for _ in range(T)]
    if kind == "random":
        return [(rng.randint(0, MASK64), rng.randint(0, MASK64)) for _ in range(T)]
    if kind == "edge_e":
        # e = 0, 1, 2, 2^k, 2^64-1 など特殊値を含む
        out = []
        for _ in range(T):
            a = rng.choice([0, 1, MASK64, rng.randint(0, MASK64)])
            e = rng.choice([0, 1, 2, MASK64, MASK64 - 1, 1 << rng.randint(0, 63)])
            out.append((a, e))
        return out
    raise ValueError(kind)


def write_case(name: str, T: int, kind: str) -> None:
    in_path = HERE / f"{name}.in"
    out_path = HERE / f"{name}.out"
    if in_path.exists() and out_path.exists():
        return
    print(f"gen {name} (T={T}, kind={kind})", file=sys.stderr)
    pairs = make_pairs(name, T, kind)
    with open(in_path, "w") as f:
        f.write(f"{T}\n")
        for a, e in pairs:
            f.write(f"{a} {e}\n")
    results = gf2_64.batch_pow(pairs)
    with open(out_path, "w") as f:
        for r in results:
            f.write(f"{r}\n")


def main() -> None:
    targets = sys.argv[1:] or list(CASES.keys())
    for name in targets:
        if name not in CASES:
            print(f"unknown: {name}", file=sys.stderr); continue
        T, kind = CASES[name]
        write_case(name, T, kind)


if __name__ == "__main__":
    main()
