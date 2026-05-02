#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""塔表現用の precomputed 定数を計算して C++ ヘッダ tower_custom_data.hpp に出力。

塔: F_{2^64} ≅ F_{2^16}[ω] / q(ω) where F_{2^16} = GF(2)[s]/(r(s))。
poly 基底 (= GF(2)[x]/(x^64+x^4+x^3+x+1)) との間の基底変換を計算。
"""
from __future__ import annotations

import sys
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gf2_64  # noqa: E402

MASK64 = gf2_64.MASK64
gf_mul = gf2_64.gf_mul
gf_pow = gf2_64.gf_pow

def f16_mul_with_r(a: int, b: int, r_low: int) -> int:
    p = 0
    for i in range(16):
        if (b >> i) & 1:
            p ^= a << i
    for i in range(30, 15, -1):
        if (p >> i) & 1:
            p ^= 1 << i
            p ^= r_low << (i - 16)
    return p & 0xFFFF


def f16_pow_with_r(a: int, k: int, r_low: int) -> int:
    res = 1
    while k:
        if k & 1:
            res = f16_mul_with_r(res, a, r_low)
        a = f16_mul_with_r(a, a, r_low)
        k >>= 1
    return res


def is_irreducible(r_low: int) -> bool:
    """r(s) が GF(2)[s] 上で既約か簡易チェック。
    deg-16 では x^(2^k) - x の各 k=1..8 について gcd を取るのが正攻法だが、
    ここでは "x が deg 16 over GF(2)" のチェックで近似。"""
    # x^(2^16) ≡ x (mod r) は常に成り立つ。さらに x^(2^k) ≠ x for k=1..8 (= proper divisors check)
    # 16 の真の約数 = 1, 2, 4, 8 → x^(2^k) ≠ x をそれぞれ確認すれば既約
    x = 2
    for k in (1, 2, 4, 8):
        # x^(2^k) を計算 (mod r)
        v = x
        for _ in range(k):
            v = f16_mul_with_r(v, v, r_low)
        if v == x:
            return False
    return True


def is_primitive_r(r_low: int) -> bool:
    """r が primitive (= x が F_{2^16}^* の生成元) か。"""
    if not is_irreducible(r_low):
        return False
    # 65535 = 3 * 5 * 17 * 257
    for p in (3, 5, 17, 257):
        if f16_pow_with_r(2, 65535 // p, r_low) == 1:
            return False
    return True


# r(s) = s^16 + s^12 + s^3 + s + 1 (CRC-16-CCITT 互換、よく使われる primitive 既約)。
R_LOW = (1 << 12) | (1 << 3) | (1 << 1) | 1  # = 0x100B
assert is_primitive_r(R_LOW), "r is not primitive"
print(f"r(s) = s^16 + 0x{R_LOW:04x} (low part)", file=sys.stderr)


def f16_mul(a: int, b: int) -> int:
    return f16_mul_with_r(a, b, R_LOW)


def f16_pow(a: int, k: int) -> int:
    return f16_pow_with_r(a, k, R_LOW)


def eval_r_at(x: int) -> int:
    """r(x) を F_{2^64} (poly 基底) で評価。"""
    val = gf_pow(x, 16)
    for i in range(16):
        if (R_LOW >> i) & 1:
            val ^= gf_pow(x, i)
    return val


# σ ∈ F_{2^64} with r(σ) = 0 を見つける。
# σ は order 65535 を持つ。F_{2^64}^* で order 65535 の元は g^((2^64-1)/65535) の形 (g は primitive)。
# poly 基底で g = 2 が primitive (= P が primitive 既約)。
G64 = 2
SIGMA_BASE = gf_pow(G64, MASK64 // 65535)

SIGMA = None
for j in range(1, 65535):
    if gcd(j, 65535) != 1:
        continue
    cand = gf_pow(SIGMA_BASE, j)
    if eval_r_at(cand) == 0:
        SIGMA = cand
        break
assert SIGMA is not None, "no σ found"
print(f"σ = 0x{SIGMA:016x}", file=sys.stderr)

# ω = G64 = 2。これは F_{2^64}^* の primitive なので degree 4 over F_{2^16}。
OMEGA = G64

# 基底 (16i + j → σ^j ω^i 互換にして u64 の bit 配置 a_i が ω^i の係数になるよう)
# bit_index k = 16*i + j (0 ≤ i ≤ 3, 0 ≤ j ≤ 15) → σ^j ω^i
basis = [0] * 64
for i in range(4):
    omega_i = gf_pow(OMEGA, i)
    for j in range(16):
        sigma_j = gf_pow(SIGMA, j)
        basis[16 * i + j] = gf_mul(sigma_j, omega_i)

# 64x64 行列 M: 列 col = basis[col] (poly 基底 u64) の bit 列
# M[row][col] = (basis[col] >> row) & 1
M = [[0] * 64 for _ in range(64)]
for col in range(64):
    for row in range(64):
        M[row][col] = (basis[col] >> row) & 1


def invert_gf2(A: list[list[int]]) -> list[list[int]]:
    n = len(A)
    M = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(A)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if M[r][col]), None)
        assert pivot is not None
        M[col], M[pivot] = M[pivot], M[col]
        for r in range(n):
            if r != col and M[r][col]:
                for c in range(2 * n):
                    M[r][c] ^= M[col][c]
    return [row[n:] for row in M]


M_inv = invert_gf2(M)
print("computed M^{-1}", file=sys.stderr)


def byte_tables(M_64x64: list[list[int]]) -> list[list[int]]:
    """各 byte 位置 (0..7) について、入力 byte 値 (0..255) → 出力 u64 のテーブル。"""
    tables = []
    for byte_pos in range(8):
        tbl = [0] * 256
        for byte_val in range(256):
            v = 0
            for bit in range(8):
                if (byte_val >> bit) & 1:
                    col = byte_pos * 8 + bit
                    for row in range(64):
                        if M_64x64[row][col]:
                            v ^= 1 << row
            tbl[byte_val] = v
        tables.append(tbl)
    return tables


# tower → poly: 入力 (tower 基底 u64) の bit i → poly 基底 vector (basis[i]) を XOR
# = M を適用する (M[row][col] = basis vector の col 番目の bit row)
TOWER_TO_POLY_BYTE = byte_tables(M)
# poly → tower: M^{-1} を適用
POLY_TO_TOWER_BYTE = byte_tables(M_inv)
print("computed byte tables", file=sys.stderr)


def apply_byte_table(tables: list[list[int]], x: int) -> int:
    """GF(2) ベクトルを byte table 経由で線形変換 (XOR で集約)。"""
    out = 0
    for pos in range(8):
        out ^= tables[pos][(x >> (pos * 8)) & 0xFF]
    return out


# q の係数。q(y) = (y + ω)(y + ω^{2^16})(y + ω^{2^32})(y + ω^{2^48}) を展開。
# Galois 共役の積は F_{2^16} に落ちる。
roots = [gf_pow(OMEGA, 1 << k) for k in (0, 16, 32, 48)]


def expand_product(rs: list[int]) -> list[int]:
    """∏(y + r_k) を多項式 [c_0, c_1, ..., c_n] (低次から) で返す。係数は F_{2^64} 上。"""
    p = [1]  # 1 (deg 0)
    for r in rs:
        # 新しい p' = (y + r) * p = shift p で y 倍 + r * p の XOR
        new = [0] * (len(p) + 1)
        for i, c in enumerate(p):
            new[i + 1] ^= c
            new[i] ^= gf_mul(r, c)
        p = new
    return p


q_poly = expand_product(roots)  # length 5: y^4 + a_3 y^3 + a_2 y^2 + a_1 y + a_0
assert len(q_poly) == 5 and q_poly[4] == 1, "q is not monic deg 4"
# a_i は F_{2^16} ⊂ F_{2^64} (poly 基底) → tower 基底に変換すると下位 16 bit に落ちる
Q_COEF_F16 = []
for i in range(4):
    a_i_poly = q_poly[i]
    a_i_tower = apply_byte_table(POLY_TO_TOWER_BYTE, a_i_poly)
    assert a_i_tower >> 16 == 0, f"a_{i} not in F_{{2^16}}: tower = 0x{a_i_tower:016x}"
    Q_COEF_F16.append(a_i_tower & 0xFFFF)
print(f"q coefs: a0={Q_COEF_F16[0]:#x}, a1={Q_COEF_F16[1]:#x}, "
      f"a2={Q_COEF_F16[2]:#x}, a3={Q_COEF_F16[3]:#x}", file=sys.stderr)


# 出力
def emit_byte_table(name: str, tbl: list[list[int]]) -> str:
    lines = [f"alignas(64) inline constexpr uint64_t {name}[8][256] = {{"]
    for pos, row in enumerate(tbl):
        lines.append(f"  {{")
        for chunk_start in range(0, 256, 4):
            chunk = ", ".join(f"0x{row[i]:016x}ull" for i in range(chunk_start, chunk_start + 4))
            lines.append(f"    {chunk},")
        lines.append(f"  }},")
    lines.append("};")
    return "\n".join(lines)


out = HERE / "tower_custom_data.hpp"
with open(out, "w") as f:
    f.write("// Auto-generated by gen_tower_custom.py. DO NOT EDIT.\n")
    f.write("#pragma once\n#include <cstdint>\n\n")
    f.write("namespace gf2_64_tower_custom {\n\n")
    f.write(f"// r(s) = s^16 + s^5 + s^3 + s + 1 (primitive irreducible)\n")
    f.write(f"inline constexpr uint16_t R_LOW = 0x{R_LOW:04x};\n\n")
    f.write(f"// σ ∈ F_{{2^64}} (poly 基底) with r(σ) = 0 (= F_{{2^16}}^* の generator)\n")
    f.write(f"inline constexpr uint64_t SIGMA_POLY = 0x{SIGMA:016x}ull;\n")
    f.write(f"inline constexpr uint64_t OMEGA_POLY = 0x{OMEGA:016x}ull;\n\n")
    f.write(f"// q(y) = y^4 + Q_COEF[3] y^3 + Q_COEF[2] y^2 + Q_COEF[1] y + Q_COEF[0]\n")
    f.write(f"inline constexpr uint16_t Q_COEF[4] = {{0x{Q_COEF_F16[0]:04x}, 0x{Q_COEF_F16[1]:04x}, 0x{Q_COEF_F16[2]:04x}, 0x{Q_COEF_F16[3]:04x}}};\n\n")
    f.write("// poly 基底 → tower 基底 byte tables\n")
    f.write(emit_byte_table("POLY_TO_TOWER_BYTE", POLY_TO_TOWER_BYTE) + "\n\n")
    f.write("// tower 基底 → poly 基底 byte tables\n")
    f.write(emit_byte_table("TOWER_TO_POLY_BYTE", TOWER_TO_POLY_BYTE) + "\n\n")
    f.write("} // namespace gf2_64_tower_custom\n")
print(f"wrote {out}", file=sys.stderr)
