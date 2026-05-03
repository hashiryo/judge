#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""yosupo discrete_logarithm_fixed_mod 用テストケース生成。
p = 999999937 (固定の素数)、g (primitive root) を選び、random query を出力。
答え (.out) は naive で BSGS して与える。
"""
from __future__ import annotations
import random
import sys
from pathlib import Path

P = 999999937
PHI = P - 1

def is_primitive_root(g: int, p: int, phi: int, phi_factors: list[int]) -> bool:
    for q in phi_factors:
        if pow(g, phi // q, p) == 1:
            return False
    return True

# phi(p) = p - 1 = 999999936 = 2^6 * 3 * 13 * 31 * 263 * 311 (たぶん)
# ファクタ計算
def factor(n: int) -> list[int]:
    fs = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            fs.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: fs.append(n)
    return fs

phi_factors = factor(PHI)

# 適当な primitive root を探す
g = 2
while not is_primitive_root(g, P, PHI, phi_factors):
    g += 1

random.seed(42)
N = 1000  # T=1000 (random_00)
queries = [random.randint(0, P - 1) for _ in range(N)]

HERE = Path(__file__).resolve().parent
inp = HERE / "random_00.in"
out = HERE / "random_00.out"

with open(inp, "w") as f:
    f.write(f"{P} {g} {N}\n")
    for q in queries: f.write(f"{q}\n")

# 答えは naive な pow で計算 (= log_g(q) を BSGS で出すのは Python では遅いので、
#   テスト用に generate-and-verify アプローチ: 各 q について x をランダム選び q = g^x mod p、
#   x を答えとする)
# だが random query なので x が未知。BSGS or 専用算法でないと求まらない。
# シンプルな解決: random query を生成する代わりに、random x から q = g^x を作る形に変更。
queries2 = []
answers = []
for _ in range(N):
    x = random.randint(0, PHI - 1)
    q = pow(g, x, P)
    # ランダムに 0 を混ぜる (~5% の確率で)
    if random.random() < 0.05:
        q = 0
        ans = -1
    else:
        ans = x
    queries2.append(q)
    answers.append(ans)

# 上書き
with open(inp, "w") as f:
    f.write(f"{P} {g} {N}\n")
    for q in queries2: f.write(f"{q}\n")
with open(out, "w") as f:
    for a in answers: f.write(f"{a}\n")

print(f"Generated {N} queries, p={P}, g={g}", file=sys.stderr)
