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

// tower mul: F_{2^16}[ω]/q(ω) 上の積
// (a_3 ω^3 + a_2 ω^2 + a_1 ω + a_0)(b_3 ω^3 + b_2 ω^2 + b_1 ω + b_0)
// = c_6 ω^6 + c_5 ω^5 + c_4 ω^4 + c_3 ω^3 + c_2 ω^2 + c_1 ω + c_0
// reduce: ω^4 = Q_COEF[3] ω^3 + Q_COEF[2] ω^2 + Q_COEF[1] ω + Q_COEF[0]
//         ω^5 = ω · ω^4
//         ω^6 = ω · ω^5
[[gnu::always_inline]] inline u64 tower_mul(u64 A, u64 B) {
 const u16 a0 = u16(A), a1 = u16(A >> 16), a2 = u16(A >> 32), a3 = u16(A >> 48);
 const u16 b0 = u16(B), b1 = u16(B >> 16), b2 = u16(B >> 32), b3 = u16(B >> 48);
 // schoolbook 16 mul
 u16 c0 = f16_mul(a0, b0);
 u16 c1 = f16_mul(a0, b1) ^ f16_mul(a1, b0);
 u16 c2 = f16_mul(a0, b2) ^ f16_mul(a1, b1) ^ f16_mul(a2, b0);
 u16 c3 = f16_mul(a0, b3) ^ f16_mul(a1, b2) ^ f16_mul(a2, b1) ^ f16_mul(a3, b0);
 u16 c4 = f16_mul(a1, b3) ^ f16_mul(a2, b2) ^ f16_mul(a3, b1);
 u16 c5 = f16_mul(a2, b3) ^ f16_mul(a3, b2);
 u16 c6 = f16_mul(a3, b3);

 // reduce (ω^4, ω^5, ω^6) を (ω^0..ω^3) に展開して c0..c3 に xor
 const u16 q0 = Q_COEF[0], q1 = Q_COEF[1], q2 = Q_COEF[2], q3 = Q_COEF[3];
 // ω^4 = q3 ω^3 + q2 ω^2 + q1 ω + q0
 // ω^5 = ω · ω^4 = q0 ω + q1 ω^2 + q2 ω^3 + q3 ω^4
 //     = q0 ω + q1 ω^2 + q2 ω^3 + q3 (q0 + q1 ω + q2 ω^2 + q3 ω^3)
 //     = (q3 q0) + (q0 + q3 q1) ω + (q1 + q3 q2) ω^2 + (q2 + q3 q3) ω^3
 // ω^6 = ω · ω^5 = (q3 q0) ω + (q0 + q3 q1) ω^2 + (q1 + q3 q2) ω^3 + (q2 + q3 q3) ω^4
 //     = (q3 q0) ω + (q0 + q3 q1) ω^2 + (q1 + q3 q2) ω^3 + (q2 + q3 q3)(q0 + q1 ω + q2 ω^2 + q3 ω^3)
 //     = (q2 + q3^2) q0 + ((q3 q0) + (q2 + q3^2) q1) ω + ...
 //
 // 直接書くより、c6 → c5+ω·c6 (= c5'); c5' → c4+ω·c5' (= c4'); c4' → c0..c3 補正、で 1 段ずつ降ろす方が速い。
 // 但しここでは "ω^4 = q0+q1 ω+q2 ω^2+q3 ω^3" の関係を使う:
 //   ω · v_3 ω^3 = v_3 q0 + v_3 q1 ω + v_3 q2 ω^2 + v_3 q3 ω^3
 // つまり「ω 倍 = shift by 1 + ω^4 への補正」
 auto omega_times = [&](u16 v0, u16 v1, u16 v2, u16 v3) {
  // 入力: v_0 + v_1 ω + v_2 ω^2 + v_3 ω^3 (∈ F_{2^16}[ω], deg < 4)
  // 出力: ω · 入力 = v_0 ω + v_1 ω^2 + v_2 ω^3 + v_3 ω^4
  //              = (v_3 q0) + (v_0 + v_3 q1) ω + (v_1 + v_3 q2) ω^2 + (v_2 + v_3 q3) ω^3
  return std::array<u16, 4>{
   f16_mul(v3, q0),
   u16(v0 ^ f16_mul(v3, q1)),
   u16(v1 ^ f16_mul(v3, q2)),
   u16(v2 ^ f16_mul(v3, q3)),
  };
 };

 // c4..c6 を順に c0..c3 に降ろす
 // c6 ω^6 = ω^2 · c6 ω^4 = ω^2 · (c6 q0 + c6 q1 ω + c6 q2 ω^2 + c6 q3 ω^3)
 // 段階的に: c6_omega = ω · (c6 q0, c6 q1, c6 q2, c6 q3)
 //           c6_omega2 = ω · c6_omega
 //           c0..c3 ^= c6_omega2
 // 同様に c5 ω^5 = ω · c5 ω^4 → ω 倍 1 回
 //         c4 ω^4 = c4 ω^4 → 0 回 (直接補正)
 std::array<u16, 4> r = {f16_mul(c4, q0), u16(f16_mul(c4, q1)), u16(f16_mul(c4, q2)), u16(f16_mul(c4, q3))};
 // c5: ω · (c5 q0, c5 q1, c5 q2, c5 q3) を加算
 {
  std::array<u16, 4> v = {f16_mul(c5, q0), f16_mul(c5, q1), f16_mul(c5, q2), f16_mul(c5, q3)};
  v = omega_times(v[0], v[1], v[2], v[3]);
  r[0] ^= v[0]; r[1] ^= v[1]; r[2] ^= v[2]; r[3] ^= v[3];
 }
 // c6: ω^2 · (c6 q0, c6 q1, c6 q2, c6 q3) を加算
 {
  std::array<u16, 4> v = {f16_mul(c6, q0), f16_mul(c6, q1), f16_mul(c6, q2), f16_mul(c6, q3)};
  v = omega_times(v[0], v[1], v[2], v[3]);
  v = omega_times(v[0], v[1], v[2], v[3]);
  r[0] ^= v[0]; r[1] ^= v[1]; r[2] ^= v[2]; r[3] ^= v[3];
 }

 c0 ^= r[0]; c1 ^= r[1]; c2 ^= r[2]; c3 ^= r[3];

 return u64(c0) | (u64(c1) << 16) | (u64(c2) << 32) | (u64(c3) << 48);
}

[[gnu::always_inline]] inline u64 mul_via_tower(u64 a_poly, u64 b_poly) {
 u64 a_t = poly_to_tower(a_poly);
 u64 b_t = poly_to_tower(b_poly);
 u64 c_t = tower_mul(a_t, b_t);
 return tower_to_poly(c_t);
}

} // namespace gf2_64_tower_custom
