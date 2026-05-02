#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""塔表現 (sparse 版): ω を探索して q をできるだけ疎にする。

tower_custom.py との違い:
  - ω を g^k で k=1..50000 試し、min poly の非零係数数を最小化
  - 出力先は tower_custom_sparse_data.hpp、namespace は gf2_64_tower_custom_sparse_sparse

通常 q は 5-term (全 4 係数 nonzero) だが、運が良ければ 4-term や 3-term q が
得られて reduce のコストが減る (1-2 mul 削減)。
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


def f16_mul_via_table(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return _F16_PW[(_F16_LN[a] + _F16_LN[b]) % 65535]


def f16_pow_via_table(a: int, k: int) -> int:
    if a == 0:
        return 0 if k > 0 else 1
    return _F16_PW[(_F16_LN[a] * k) % 65535]


# log/exp テーブルを Python 側でも 1 度だけ構築 (検索高速化)
_F16_PW = [0] * 65536
_F16_LN = [0] * 65536
_cur = 1
for _k in range(65535):
    _F16_PW[_k] = _cur
    _F16_LN[_cur] = _k
    _cur = f16_mul_with_r(_cur, 2, R_LOW)
_F16_PW[65535] = 1


def f16_mul(a: int, b: int) -> int:
    return f16_mul_via_table(a, b)


def f16_pow(a: int, k: int) -> int:
    return f16_pow_via_table(a, k)


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


# =============================================================================
# 3-項 q の探索 + 対応する ω を F_2 線形方程式で求める
# =============================================================================
def f16_to_f64(v: int) -> int:
    """F_{2^16} 元 v (= s-基底の 16-bit) を F_{2^64} の poly 基底に埋め込む。
    s = σ という対応で、bit i が立つなら σ^i を加算。"""
    out = 0
    for i in range(16):
        if (v >> i) & 1:
            out ^= gf_pow(SIGMA, i)
    return out


def f16_poly_squaring_mod_q(cur: list[int], q: list[int]) -> list[int]:
    """F_{2^16}[y] / q(y) で cur^2 を計算 (deg < 4 を保つ)。
    cur, q は係数リスト (低次から)。q は monic deg 4。"""
    # squaring: char 2 で (sum c_i y^i)^2 = sum c_i^2 y^{2i}
    new = [0] * (2 * len(cur) - 1)
    for i, c in enumerate(cur):
        new[2 * i] = f16_mul(c, c)
    # 末尾を高次から消去 (q[4] = 1 monic 前提で y^d → -q の rewrite)
    for d in range(len(new) - 1, 3, -1):
        if new[d] == 0:
            continue
        top = new[d]
        for j in range(5):
            new[(d - 4) + j] ^= f16_mul(top, q[j])
    return new[:4]


def has_f16_root(beta_f16: int, alpha_f16: int) -> bool:
    """y^4 + β y + α が F_{2^16} に根を持つか (brute force)。"""
    for y in range(65536):
        val = f16_mul(f16_mul(y, y), f16_mul(y, y))  # y^4
        val ^= f16_mul(beta_f16, y)
        val ^= alpha_f16
        if val == 0:
            return True
    return False


def f16_inv(a: int) -> int:
    if a == 0:
        raise ZeroDivisionError("0 inv")
    return _F16_PW[(65535 - _F16_LN[a]) % 65535]


def f16_poly_mod(a: list[int], b: list[int]) -> list[int]:
    """a mod b in F_{2^16}[y]. b の最高次係数 nonzero 前提。"""
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    deg_b = len(b) - 1
    while len(a) > deg_b:
        if a[-1] == 0:
            a.pop()
            continue
        scale = f16_mul(a[-1], f16_inv(b[-1]))
        shift = len(a) - 1 - deg_b
        for i, c in enumerate(b):
            a[i + shift] ^= f16_mul(scale, c)
        while a and a[-1] == 0:
            a.pop()
    return a


def f16_poly_gcd(a: list[int], b: list[int]) -> list[int]:
    a = [c for c in a]
    b = [c for c in b]
    while a and a[-1] == 0:
        a.pop()
    while b and b[-1] == 0:
        b.pop()
    while b:
        a, b = b, f16_poly_mod(a, b)
        while b and b[-1] == 0:
            b.pop()
    return a


def is_q_irreducible_over_f16(q: list[int]) -> bool:
    """q (monic deg 4 in F_{2^16}[y]) が既約か。
    gcd(q, y^{2^16} - y) = 1 (F_{2^16}-root 無し) かつ
    gcd(q, y^{2^32} - y) = 1 (F_{2^32}-root 無し = 二次因子無し) を確認。
    両方満たせば deg-4 で既約 (cubic 因子は F_{2^48} 系で F_{2^64} に現れない)。"""
    # y^{2^16} mod q
    cur = [0, 1, 0, 0]
    for _ in range(16):
        cur = f16_poly_squaring_mod_q(cur, q)
    # diff = cur + y (char 2)
    diff = list(cur)
    if len(diff) <= 1:
        diff = (diff + [0] * 2)[:2]
    diff[1] ^= 1
    while diff and diff[-1] == 0:
        diff.pop()
    if not diff:
        return False  # y^{2^16} = y
    g = f16_poly_gcd(q, diff)
    if len(g) > 1:
        return False  # gcd has degree ≥ 1 → F_{2^16}-root が存在
    # y^{2^32} mod q
    cur2 = list(cur)
    for _ in range(16):
        cur2 = f16_poly_squaring_mod_q(cur2, q)
    diff2 = list(cur2)
    if len(diff2) <= 1:
        diff2 = (diff2 + [0] * 2)[:2]
    diff2[1] ^= 1
    while diff2 and diff2[-1] == 0:
        diff2.pop()
    if not diff2:
        return False
    g2 = f16_poly_gcd(q, diff2)
    if len(g2) > 1:
        return False  # 二次因子が存在
    return True


def find_omega_root(beta_f16: int, alpha_f16: int) -> int | None:
    """ω^4 + β·ω + α = 0 を F_{2^64} の poly 基底で解く。
    ω → ω^4 と ω → β·ω は GF(2)-線形なので、64×64 線形方程式。"""
    beta_f64 = f16_to_f64(beta_f16)
    alpha_f64 = f16_to_f64(alpha_f16)
    # M = ω → ω^4 + β·ω の行列表示 (64x64 over F_2)
    M = [[0] * 64 for _ in range(64)]
    for col in range(64):
        e = 1 << col
        e2 = gf_mul(e, e)
        e4 = gf_mul(e2, e2)
        be = gf_mul(beta_f64, e)
        result = e4 ^ be
        for row in range(64):
            M[row][col] = (result >> row) & 1
    # 拡大係数行列に α_f64 を追加して Gauss 消去
    aug = [M[i] + [(alpha_f64 >> i) & 1] for i in range(64)]
    pivot_cols = []
    row = 0
    for col in range(64):
        # ピボット候補を探す
        p = None
        for r in range(row, 64):
            if aug[r][col]:
                p = r
                break
        if p is None:
            continue
        aug[row], aug[p] = aug[p], aug[row]
        for r in range(64):
            if r != row and aug[r][col]:
                for c in range(65):
                    aug[r][c] ^= aug[row][c]
        pivot_cols.append(col)
        row += 1
    # 解の存在確認: ピボット行のない行で右辺が立っていれば解なし
    for r in range(row, 64):
        if aug[r][64] != 0:
            return None
    # 自由変数を 0 に設定して 1 つの解を取り出す
    omega_bits = [0] * 64
    for i, p_col in enumerate(pivot_cols):
        omega_bits[p_col] = aug[i][64]
    return sum(b << i for i, b in enumerate(omega_bits))


# 3-項 q = y^4 + β y + α の探索 (β, α ∈ F_{2^16}, α ≠ 0)。
# 既約な q を見つける確率が低そう (= 順次探索が遅い) ので、まず random で試す。
import random as _random
_random.seed(42)

def expand_product_f64(roots: list[int]) -> list[int]:
    """∏(y + r_k) を多項式係数 (低次から) で返す。係数は F_{2^64} 上。"""
    p = [1]
    for r in roots:
        new = [0] * (len(p) + 1)
        for i, c in enumerate(p):
            new[i + 1] ^= c
            new[i] ^= gf_mul(r, c)
        p = new
    return p


def min_poly_of(omega_val: int) -> list[int]:
    """ω の F_{2^16} 上の minimal poly (係数は F_{2^64} の poly 基底 u64 で表現)。"""
    rs = [gf_pow(omega_val, 1 << (16 * i)) for i in range(4)]
    return expand_product_f64(rs)


print("searching ω with sparse min poly (this may take a few minutes)...", file=sys.stderr)
best_nonzero = 5
best_omega = G64  # fallback ω = 2
for k_test in range(1, 50000):
    omega_test = gf_pow(G64, k_test)
    # degree 4 over F_{2^16} を確認 (= F_{2^32} に落ちない)
    if gf_pow(omega_test, 1 << 32) == omega_test:
        continue
    q_f64_test = min_poly_of(omega_test)
    nonzero = sum(1 for c in q_f64_test[:4] if c != 0)
    if nonzero < best_nonzero:
        best_nonzero = nonzero
        best_omega = omega_test
        print(f"  k={k_test}: nonzero={nonzero}, ω=0x{omega_test:016x}", file=sys.stderr)
        if nonzero <= 2:
            break
OMEGA = best_omega
print(f"selected ω = 0x{OMEGA:016x} (nonzero coef count = {best_nonzero})", file=sys.stderr)

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


# ω^4, ω^5, ω^6 を (1, ω, ω^2, ω^3) の F_{2^16} 線形結合で表す係数を事前計算。
# ω^4 = a0 + a1 ω + a2 ω^2 + a3 ω^3 (= Q_COEF)
# ω^5 = ω · ω^4
# ω^6 = ω · ω^5
def omega_times(v: list[int]) -> list[int]:
    """ω · (v0 + v1 ω + v2 ω^2 + v3 ω^3) を返す。
    = (v3 a0) + (v0 + v3 a1) ω + (v1 + v3 a2) ω^2 + (v2 + v3 a3) ω^3"""
    a0, a1, a2, a3 = Q_COEF_F16
    return [
        f16_mul(v[3], a0),
        v[0] ^ f16_mul(v[3], a1),
        v[1] ^ f16_mul(v[3], a2),
        v[2] ^ f16_mul(v[3], a3),
    ]


OMEGA4_COEF = list(Q_COEF_F16)              # ω^4
OMEGA5_COEF = omega_times(OMEGA4_COEF)      # ω^5
OMEGA6_COEF = omega_times(OMEGA5_COEF)      # ω^6
print(f"ω^4 coefs: {[hex(c) for c in OMEGA4_COEF]}", file=sys.stderr)
print(f"ω^5 coefs: {[hex(c) for c in OMEGA5_COEF]}", file=sys.stderr)
print(f"ω^6 coefs: {[hex(c) for c in OMEGA6_COEF]}", file=sys.stderr)


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


out = HERE / "tower_custom_sparse_data.hpp"
with open(out, "w") as f:
    f.write("// Auto-generated by gen_tower_custom.py. DO NOT EDIT.\n")
    f.write("#pragma once\n#include <cstdint>\n\n")
    f.write("namespace gf2_64_tower_custom_sparse {\n\n")
    f.write(f"// r(s) = s^16 + s^5 + s^3 + s + 1 (primitive irreducible)\n")
    f.write(f"inline constexpr uint16_t R_LOW = 0x{R_LOW:04x};\n\n")
    f.write(f"// σ ∈ F_{{2^64}} (poly 基底) with r(σ) = 0 (= F_{{2^16}}^* の generator)\n")
    f.write(f"inline constexpr uint64_t SIGMA_POLY = 0x{SIGMA:016x}ull;\n")
    f.write(f"inline constexpr uint64_t OMEGA_POLY = 0x{OMEGA:016x}ull;\n\n")
    f.write(f"// q(y) = y^4 + Q_COEF[3] y^3 + Q_COEF[2] y^2 + Q_COEF[1] y + Q_COEF[0]\n")
    f.write(f"inline constexpr uint16_t Q_COEF[4] = {{0x{Q_COEF_F16[0]:04x}, 0x{Q_COEF_F16[1]:04x}, 0x{Q_COEF_F16[2]:04x}, 0x{Q_COEF_F16[3]:04x}}};\n\n")
    f.write(f"// ω^k を (1, ω, ω^2, ω^3) の F_{{2^16}} 線形結合で表した係数 (k=4,5,6)。\n")
    f.write(f"// reduce で c_k × OMEGA{{k}}_COEF[i] を c0..c3 に xor すれば良い。\n")
    f.write(f"inline constexpr uint16_t OMEGA4_COEF[4] = {{0x{OMEGA4_COEF[0]:04x}, 0x{OMEGA4_COEF[1]:04x}, 0x{OMEGA4_COEF[2]:04x}, 0x{OMEGA4_COEF[3]:04x}}};\n")
    f.write(f"inline constexpr uint16_t OMEGA5_COEF[4] = {{0x{OMEGA5_COEF[0]:04x}, 0x{OMEGA5_COEF[1]:04x}, 0x{OMEGA5_COEF[2]:04x}, 0x{OMEGA5_COEF[3]:04x}}};\n")
    f.write(f"inline constexpr uint16_t OMEGA6_COEF[4] = {{0x{OMEGA6_COEF[0]:04x}, 0x{OMEGA6_COEF[1]:04x}, 0x{OMEGA6_COEF[2]:04x}, 0x{OMEGA6_COEF[3]:04x}}};\n\n")
    f.write("// poly 基底 → tower 基底 byte tables\n")
    f.write(emit_byte_table("POLY_TO_TOWER_BYTE", POLY_TO_TOWER_BYTE) + "\n\n")
    f.write("// tower 基底 → poly 基底 byte tables\n")
    f.write(emit_byte_table("TOWER_TO_POLY_BYTE", TOWER_TO_POLY_BYTE) + "\n\n")
    f.write("} // namespace gf2_64_tower_custom_sparse\n")
print(f"wrote {out}", file=sys.stderr)
