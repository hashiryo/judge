#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""GF(2^64) 離散対数 log_g(x) のテストケース生成。

g = 2 は P(x)=x^64+x^4+x^3+x+1 の下で原始元 (位数 2^64-1)。
入力 x は g^k (k は random in [0, 2^64-2]) で構成。x ≠ 0。
期待出力 (.out) は k だが、checker.cpp が g^k_user == x で再検証するので
複数解可。
"""
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "_shared"))
import gf2_64  # noqa: E402

MASK64 = gf2_64.MASK64
ORDER = gf2_64.ORDER  # = 2^64 - 1
G = gf2_64.PRIMITIVE_G  # = 2

CASES = {
    "sample_00":           (8,    "sample"),
    "small_00":            (50,   "small"),
    "edge_00":             (50,   "edge"),
    # log は BSGS で 1 query ~数 ms 〜数十 ms かかるので T は控えめに
    "random_medium_00":    (1000, "random"),
    "random_00":           (10000, "random"),
}


def make_ks(name: str, T: int, kind: str) -> list[int]:
    """正解 k のリストを返す。x は後で batch_log_inv (g^k) で一括計算。"""
    rng = random.Random(hash(("log", name)) & 0xFFFFFFFF)
    if kind == "sample":
        return [0, 1, 2, 100, 1 << 16, 1 << 32, 1 << 60, ORDER - 1][:T]
    if kind == "small":
        return [rng.randint(0, 1000) for _ in range(T)]
    if kind == "edge":
        return [
            rng.choice([0, 1, 2, ORDER - 2, ORDER - 1, rng.randint(0, ORDER - 1)])
            for _ in range(T)
        ]
    if kind == "random":
        return [rng.randint(0, ORDER - 1) for _ in range(T)]
    raise ValueError(kind)


def write_case(name: str, T: int, kind: str) -> None:
    in_path = HERE / f"{name}.in"
    out_path = HERE / f"{name}.out"
    if in_path.exists() and out_path.exists():
        return
    print(f"gen {name} (T={T}, kind={kind})", file=sys.stderr)
    ks = make_ks(name, T, kind)
    xs = gf2_64.batch_log_inv(ks)  # x_i = g^{k_i}
    pairs = list(zip(xs, ks))
    with open(in_path, "w") as f:
        f.write(f"{T}\n")
        for x, _k in pairs:
            f.write(f"{x}\n")
    with open(out_path, "w") as f:
        for _x, k in pairs:
            f.write(f"{k}\n")


def main() -> None:
    targets = sys.argv[1:] or list(CASES.keys())
    for name in targets:
        if name not in CASES:
            print(f"unknown: {name}", file=sys.stderr); continue
        T, kind = CASES[name]
        write_case(name, T, kind)


if __name__ == "__main__":
    main()
