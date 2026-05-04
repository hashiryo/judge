#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""GF(2^64) で a^4 を計算するテストケース生成。"""
import random
import sys
from pathlib import Path

HERE= Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "_shared"))
import gf2_64  # noqa: E402

MASK64= gf2_64.MASK64

CASES: dict[str, tuple[int, str]]= {
    "sample_00":           (8,    "sample"),
    "small_00":            (100,  "small"),
    "edge_zero_one_00":    (200,  "edge_zero_one"),
    "structured_00":       (1000, "structured"),
    "random_00":           (10000, "random"),
    "random_large_00":     (100_000, "random"),
}


def make_values(name: str, T: int, kind: str) -> list[int]:
    rng= random.Random(hash(("frob4", name)) & 0xFFFFFFFF)
    if kind == "sample":
        return [0, 1, 2, 3, 0xFF, MASK64,
                0x12345678ABCDEF00, 0xFEDCBA9876543210][:T]
    if kind == "small":
        return [rng.randint(0, 255) for _ in range(T)]
    if kind == "random":
        return [rng.randint(0, MASK64) for _ in range(T)]
    if kind == "edge_zero_one":
        return [rng.choice([0, 1, MASK64, rng.randint(0, MASK64)]) for _ in range(T)]
    if kind == "structured":
        return [(1 << rng.randint(0, 63)) | rng.randint(0, 7) for _ in range(T)]
    raise ValueError(kind)


def write_case(name: str, T: int, kind: str) -> None:
    in_path= HERE / f"{name}.in"
    out_path= HERE / f"{name}.out"
    if in_path.exists() and out_path.exists():
        return
    print(f"gen {name} (T={T}, kind={kind})", file=sys.stderr)
    values= make_values(name, T, kind)
    # a^4 = pow(a, 4) を batch_pow で
    pairs= [(a, 4) for a in values]
    results= gf2_64.batch_pow(pairs)
    with open(in_path, "w") as f:
        f.write(f"{T}\n")
        for a in values:
            f.write(f"{a}\n")
    with open(out_path, "w") as f:
        for r in results:
            f.write(f"{r}\n")


def main() -> None:
    targets= sys.argv[1:] or list(CASES.keys())
    for name in targets:
        if name not in CASES:
            print(f"unknown: {name}", file=sys.stderr)
            continue
        T, kind= CASES[name]
        write_case(name, T, kind)


if __name__ == "__main__":
    main()
