#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""mod-inv-prime のランダムテストケース生成。
p = 998244353 (NTT 素数) 固定で、N 個のクエリを random に生成。
~5% の確率で 0 を含む (sentinel テスト)。
"""
from __future__ import annotations
import random
import sys
from pathlib import Path

P = 998244353
random.seed(42)

HERE = Path(__file__).resolve().parent

def gen(name: str, N: int):
    queries = []
    for _ in range(N):
        if random.random() < 0.05:
            queries.append(0)
        else:
            queries.append(random.randint(1, P - 1))
    inp = HERE / f"{name}.in"
    out = HERE / f"{name}.out"
    with open(inp, "w") as f:
        f.write(f"{P} {N}\n")
        for q in queries: f.write(f"{q}\n")
    with open(out, "w") as f:
        for q in queries:
            if q == 0: f.write("-1\n")
            else: f.write(f"{pow(q, P - 2, P)}\n")
    print(f"Generated {name}: T={N}", file=sys.stderr)

gen("random_00", 100000)         # T=100K (init コストが効くサイズ)
gen("random_medium_00", 10000)   # T=10K  (callgrind 用)
gen("random_small_00", 100)      # T=100  (small)
