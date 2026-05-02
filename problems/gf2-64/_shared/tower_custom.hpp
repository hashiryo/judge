#pragma once
// =============================================================================
// 自前 (= 非 Nimber) 塔表現での F_{2^64} 演算。
//   F_{2^16} = GF(2)[s] / (s^16 + s^12 + s^3 + s + 1)  (primitive)
//   F_{2^64} = F_{2^16}[ω] / q(ω)  where q is min poly of ω over F_{2^16}
//   tower 基底 u64: bits 0..15 = a_0, 16..31 = a_1, 32..47 = a_2, 48..63 = a_3
//                  各 a_i ∈ F_{2^16} は GF(2)[s]/r(s) の元 (s-基底 16 bit)
// log/exp テーブルは起動時に構築 (~1 ms)。基底変換テーブルは事前計算済み定数。
// =============================================================================
#include <array>
#include <bit>
#include <cstdint>

#include "tower_custom_data.hpp"

namespace gf2_64_tower_custom {
using u64 = unsigned long long;
using u16 = unsigned short;
using u8 = unsigned char;
using u32 = unsigned;

// 起動時の F_{2^16} bit-by-bit mul (log/exp テーブル構築用のみ)
inline u16 f16_mul_slow(u16 a, u16 b) {
 u32 p = 0;
 for (int i = 0; i < 16; ++i) if ((b >> i) & 1) p ^= u32(a) << i;
 for (int i = 30; i >= 16; --i) {
  if ((p >> i) & 1) {
   p ^= u32(1) << i;
   p ^= u32(R_LOW) << (i - 16);
  }
 }
 return u16(p);
}

// log/exp テーブル (namespace level static、main で 1 度 init() を呼ぶ)
inline u16 PW[65536];
inline u16 LN[65536];

inline void init_f16_tables() {
 static bool initialized = false;
 if (initialized) return;
 initialized = true;
 // s = 2 が F_{2^16}^* の primitive (r が primitive 既約)
 u16 g = 2;
 u16 cur = 1;
 for (u32 k = 0; k < 65535; ++k) {
  PW[k] = cur;
  LN[cur] = k;
  cur = f16_mul_slow(cur, g);
 }
 PW[65535] = 1;
 LN[0] = 0;
}

[[gnu::always_inline]] inline u16 f16_mul(u16 a, u16 b) {
 if (a == 0 || b == 0) return 0;
 u32 idx = u32(LN[a]) + u32(LN[b]);
 if (idx >= 65535) idx -= 65535;
 return PW[idx];
}

// 基底変換: poly 基底 u64 ↔ tower 基底 u64
[[gnu::always_inline]] inline u64 poly_to_tower(u64 v) {
 auto vb = std::bit_cast<std::array<u8, 8>>(v);
 return POLY_TO_TOWER_BYTE[0][vb[0]] ^ POLY_TO_TOWER_BYTE[1][vb[1]]
      ^ POLY_TO_TOWER_BYTE[2][vb[2]] ^ POLY_TO_TOWER_BYTE[3][vb[3]]
      ^ POLY_TO_TOWER_BYTE[4][vb[4]] ^ POLY_TO_TOWER_BYTE[5][vb[5]]
      ^ POLY_TO_TOWER_BYTE[6][vb[6]] ^ POLY_TO_TOWER_BYTE[7][vb[7]];
}
[[gnu::always_inline]] inline u64 tower_to_poly(u64 v) {
 auto vb = std::bit_cast<std::array<u8, 8>>(v);
 return TOWER_TO_POLY_BYTE[0][vb[0]] ^ TOWER_TO_POLY_BYTE[1][vb[1]]
      ^ TOWER_TO_POLY_BYTE[2][vb[2]] ^ TOWER_TO_POLY_BYTE[3][vb[3]]
      ^ TOWER_TO_POLY_BYTE[4][vb[4]] ^ TOWER_TO_POLY_BYTE[5][vb[5]]
      ^ TOWER_TO_POLY_BYTE[6][vb[6]] ^ TOWER_TO_POLY_BYTE[7][vb[7]];
}

// ω^4..ω^6 を q で reduce し c0..c3 に降ろす共通ルーチン。
// 入出力: c[0..6] (deg 6 まで持つ多項式) を c[0..3] に縮約。
//
// 高速化: ω^k = sum_i OMEGAk_COEF[i] · ω^i (k=4,5,6) を事前計算しているので、
// reduce は c_k × OMEGAk_COEF[i] を c_i に xor していくだけ。3 × 4 = 12 muls。
// (素朴計算だと c5/c6 で ω 倍を繰り返す形で 24 muls 必要だった)
[[gnu::always_inline]] inline void reduce_omega_powers(
    u16 c0, u16 c1, u16 c2, u16 c3, u16 c4, u16 c5, u16 c6,
    u16& out0, u16& out1, u16& out2, u16& out3)
{
 out0 = c0 ^ f16_mul(c4, OMEGA4_COEF[0]) ^ f16_mul(c5, OMEGA5_COEF[0]) ^ f16_mul(c6, OMEGA6_COEF[0]);
 out1 = c1 ^ f16_mul(c4, OMEGA4_COEF[1]) ^ f16_mul(c5, OMEGA5_COEF[1]) ^ f16_mul(c6, OMEGA6_COEF[1]);
 out2 = c2 ^ f16_mul(c4, OMEGA4_COEF[2]) ^ f16_mul(c5, OMEGA5_COEF[2]) ^ f16_mul(c6, OMEGA6_COEF[2]);
 out3 = c3 ^ f16_mul(c4, OMEGA4_COEF[3]) ^ f16_mul(c5, OMEGA5_COEF[3]) ^ f16_mul(c6, OMEGA6_COEF[3]);
}

// Schoolbook: 16 個の F_{2^16} 積。
[[gnu::always_inline]] inline u64 tower_mul_schoolbook(u64 A, u64 B) {
 const u16 a0 = u16(A), a1 = u16(A >> 16), a2 = u16(A >> 32), a3 = u16(A >> 48);
 const u16 b0 = u16(B), b1 = u16(B >> 16), b2 = u16(B >> 32), b3 = u16(B >> 48);
 u16 c0 = f16_mul(a0, b0);
 u16 c1 = f16_mul(a0, b1) ^ f16_mul(a1, b0);
 u16 c2 = f16_mul(a0, b2) ^ f16_mul(a1, b1) ^ f16_mul(a2, b0);
 u16 c3 = f16_mul(a0, b3) ^ f16_mul(a1, b2) ^ f16_mul(a2, b1) ^ f16_mul(a3, b0);
 u16 c4 = f16_mul(a1, b3) ^ f16_mul(a2, b2) ^ f16_mul(a3, b1);
 u16 c5 = f16_mul(a2, b3) ^ f16_mul(a3, b2);
 u16 c6 = f16_mul(a3, b3);
 u16 r0, r1, r2, r3;
 reduce_omega_powers(c0, c1, c2, c3, c4, c5, c6, r0, r1, r2, r3);
 return u64(r0) | (u64(r1) << 16) | (u64(r2) << 32) | (u64(r3) << 48);
}

// 2-level Karatsuba: 9 個の F_{2^16} 積。
//
// a = (a0 + a1 ω) + (a2 + a3 ω) ω^2 = a_low + a_high ω^2
// b = b_low + b_high ω^2
//
// 1-level Karatsuba (deg-1 × deg-1 → deg-2):
//   (p + q ω)(r + s ω) = p r + ((p+q)(r+s) + p r + q s) ω + q s ω^2
//   = 3 個の F_{2^16} 積
//
// 2-level Karatsuba (deg-3 × deg-3 → deg-6):
//   3 つの 1-level Karatsuba で 3 × 3 = 9 個の F_{2^16} 積
[[gnu::always_inline]] inline u64 tower_mul_karatsuba(u64 A, u64 B) {
 const u16 a0 = u16(A), a1 = u16(A >> 16), a2 = u16(A >> 32), a3 = u16(A >> 48);
 const u16 b0 = u16(B), b1 = u16(B >> 16), b2 = u16(B >> 32), b3 = u16(B >> 48);

 // LL = (a0 + a1 ω)(b0 + b1 ω) = LL_0 + LL_1 ω + LL_2 ω^2
 const u16 LL_0 = f16_mul(a0, b0);
 const u16 LL_2 = f16_mul(a1, b1);
 const u16 LL_1 = u16(f16_mul(u16(a0 ^ a1), u16(b0 ^ b1)) ^ LL_0 ^ LL_2);

 // HH = (a2 + a3 ω)(b2 + b3 ω)
 const u16 HH_0 = f16_mul(a2, b2);
 const u16 HH_2 = f16_mul(a3, b3);
 const u16 HH_1 = u16(f16_mul(u16(a2 ^ a3), u16(b2 ^ b3)) ^ HH_0 ^ HH_2);

 // MM = ((a0+a2) + (a1+a3) ω)((b0+b2) + (b1+b3) ω)
 const u16 ma0 = u16(a0 ^ a2), ma1 = u16(a1 ^ a3);
 const u16 mb0 = u16(b0 ^ b2), mb1 = u16(b1 ^ b3);
 const u16 MM_0 = f16_mul(ma0, mb0);
 const u16 MM_2 = f16_mul(ma1, mb1);
 const u16 MM_1 = u16(f16_mul(u16(ma0 ^ ma1), u16(mb0 ^ mb1)) ^ MM_0 ^ MM_2);

 // Combine: c = LL + (MM + LL + HH) ω^2 + HH ω^4 (in char 2: -x = x)
 const u16 c0 = LL_0;
 const u16 c1 = LL_1;
 const u16 c2 = u16(LL_2 ^ MM_0 ^ LL_0 ^ HH_0);
 const u16 c3 = u16(MM_1 ^ LL_1 ^ HH_1);
 const u16 c4 = u16(MM_2 ^ LL_2 ^ HH_2 ^ HH_0);
 const u16 c5 = HH_1;
 const u16 c6 = HH_2;

 u16 r0, r1, r2, r3;
 reduce_omega_powers(c0, c1, c2, c3, c4, c5, c6, r0, r1, r2, r3);
 return u64(r0) | (u64(r1) << 16) | (u64(r2) << 32) | (u64(r3) << 48);
}

// デフォルトはスクールブック。Karatsuba 版を使うときは tower_mul_karatsuba を直接呼ぶ。
[[gnu::always_inline]] inline u64 tower_mul(u64 A, u64 B) { return tower_mul_schoolbook(A, B); }

[[gnu::always_inline]] inline u64 mul_via_tower(u64 a_poly, u64 b_poly) {
 u64 a_t = poly_to_tower(a_poly);
 u64 b_t = poly_to_tower(b_poly);
 u64 c_t = tower_mul(a_t, b_t);
 return tower_to_poly(c_t);
}

[[gnu::always_inline]] inline u64 mul_via_tower_karatsuba(u64 a_poly, u64 b_poly) {
 u64 a_t = poly_to_tower(a_poly);
 u64 b_t = poly_to_tower(b_poly);
 u64 c_t = tower_mul_karatsuba(a_t, b_t);
 return tower_to_poly(c_t);
}

} // namespace gf2_64_tower_custom
