#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""GF(2^64) 除算 a / b のテストケース生成。b ≠ 0。"""
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "_shared"))
import gf2_64  # noqa: E402

MASK64 = gf2_64.MASK64

CASES = {
    "sample_00":           (8,    "sample"),
    "small_00":            (100,  "small"),
    "edge_one_00":         (200,  "edge_one"),
    "structured_00":       (1000, "structured"),
    "random_00":           (10000, "random"),
    "random_large_00":     (100_000, "random"),
}


def make_pairs(name: str, T: int, kind: str) -> list[tuple[int, int]]:
    rng = random.Random(hash(("div", name)) & 0xFFFFFFFF)
    if kind == "sample":
        return [
            (0, 1), (1, 1), (2, 2), (3, 5), (5, 3),
            (MASK64, 1), (MASK64, MASK64),
            (0x12345678ABCDEF00, 0xFEDCBA9876543210),
        ][:T]

    def nz():
        v = 0
        while v == 0:
            v = rng.randint(1, MASK64)
        return v

    if kind == "small":
        return [(rng.randint(0, 255), rng.randint(1, 255)) for _ in range(T)]
    if kind == "random":
        return [(rng.randint(0, MASK64), nz()) for _ in range(T)]
    if kind == "edge_one":
        out = []
        for _ in range(T):
            a = rng.choice([0, 1, MASK64, rng.randint(0, MASK64)])
            b = rng.choice([1, MASK64, nz()])
            out.append((a, b))
        return out
    if kind == "structured":
        out = []
        for _ in range(T):
            sa, sb = rng.randint(0, 63), rng.randint(0, 63)
            a = (1 << sa) | rng.randint(0, 7)
            b = (1 << sb) | rng.randint(0, 7) | 1
            out.append((a, b))
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
        for a, b in pairs:
            f.write(f"{a} {b}\n")
    results = gf2_64.batch_div(pairs)
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
