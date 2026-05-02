#pragma once
// pclmul_norm_inv16.hpp の改善版:
// nim_to_poly は入力が「上位 48 bit = 0」と分かっているので 8 lookups のうち 2 個
// だけで済む。明示的に embed_f16 (2 lookups only) を書くことで確実に高速化。
// 同様に poly_to_nim は subfield 元なので低 16 bit しか使わない → extract_f16
// (8 lookups で 16-bit 出力) で 64-bit XOR を 16-bit XOR に圧縮。
//
// 期待: 6 lookups 節約 (per inv)、~18 cycle 改善。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"
#include "../../_shared/basis_change.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_pclmul_norm_inv16_compact {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;
using u16 = unsigned short;

[[gnu::always_inline]] inline u64 spread_bits(u32 a) {
#if defined(__BMI2__) && !defined(USE_SIMDE)
 return _pdep_u64(u64(a), 0x5555555555555555ull);
#else
 u64 x = a;
 x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x = (x | (x <<  8)) & 0x00FF00FF00FF00FFull;
 x = (x | (x <<  4)) & 0x0F0F0F0F0F0F0F0Full;
 x = (x | (x <<  2)) & 0x3333333333333333ull;
 x = (x | (x <<  1)) & 0x5555555555555555ull;
 return x;
#endif
}
PCLMUL_FN u64 sq(u64 a) {
 const u64 lo = spread_bits(u32(a));
 const u64 hi = spread_bits(u32(a >> 32));
 return reduce(__m128i{(long long) lo, (long long) hi});
}

inline u64 FROB16_BYTE[8][256];
// poly basis 64-bit (subfield 元) → 16-bit nim 表現
inline u16 EXTRACT_F16_BYTE[8][256];
// nim 16-bit → poly basis 64-bit
inline u64 EMBED_F16_BYTE[2][256];
// F_{2^16} 直接 inv テーブル (nim 基底)
inline u16 INV16[65536];
inline bool inited = false;

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 // Frobenius byte table (16 sqs = ^{2^16})
 {
  u64 col[64];
  for (int j = 0; j < 64; ++j) {
   u64 v = u64(1) << j;
   for (int k = 0; k < 16; ++k) v = sq(v);
   col[j] = v;
  }
  for (int p = 0; p < 8; ++p) {
   for (int b = 0; b < 256; ++b) {
    u64 v = 0;
    for (int bit = 0; bit < 8; ++bit) {
     if ((b >> bit) & 1) v ^= col[p * 8 + bit];
    }
    FROB16_BYTE[p][b] = v;
   }
  }
 }
 // EXTRACT_F16 = low 16 bits of BASIS_BYTE
 for (int p = 0; p < 8; ++p) {
  for (int b = 0; b < 256; ++b) {
   EXTRACT_F16_BYTE[p][b] = u16(gf2_64_basis::BASIS_BYTE[p][b]);
  }
 }
 // EMBED_F16 = INV_BYTE[0..1] (低 2 byte に対する nim_to_poly 寄与)
 for (int p = 0; p < 2; ++p) {
  for (int b = 0; b < 256; ++b) {
   EMBED_F16_BYTE[p][b] = gf2_64_basis::INV_BYTE[p][b];
  }
 }
 // F_{2^16} の log/exp を一旦構築 → INV16 を作る
 u16 PW[65536], LN[65536];
 PW[0] = PW[65535] = 1;
 for (int i = 1; i < 65535; ++i) {
  PW[i] = u16((PW[i-1] << 1) ^ (0x1681fu & u16(-(PW[i-1] >= 0x8000u))));
 }
 constexpr u16 f2n[16] = {0x0001u, 0x2827u, 0x392bu, 0x8000u, 0x20fdu, 0x4d1du, 0xde4au, 0x0a17u,
                          0x3464u, 0xe3a9u, 0x6d8du, 0x34bcu, 0xa921u, 0xa173u, 0x0ebcu, 0x0e69u};
 for (int i = 1; i < 65535; ++i) {
  u16 x = PW[i], y = 0;
  for (; x; x &= x - 1) y ^= f2n[__builtin_ctz(x)];
  PW[i] = y;
  LN[y] = u16(i);
 }
 LN[1] = 0;
 INV16[0] = 0;
 INV16[1] = 1;
 for (int v = 2; v < 65536; ++v) {
  INV16[v] = PW[(65535u - u32(LN[v])) % 65535u];
 }
}

[[gnu::always_inline]] inline u64 frob16(u64 a) {
 return FROB16_BYTE[0][u8(a)]       ^ FROB16_BYTE[1][u8(a >>  8)]
      ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)]
      ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)]
      ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}
[[gnu::always_inline]] inline u16 extract_f16(u64 N) {
 return EXTRACT_F16_BYTE[0][u8(N)]       ^ EXTRACT_F16_BYTE[1][u8(N >>  8)]
      ^ EXTRACT_F16_BYTE[2][u8(N >> 16)] ^ EXTRACT_F16_BYTE[3][u8(N >> 24)]
      ^ EXTRACT_F16_BYTE[4][u8(N >> 32)] ^ EXTRACT_F16_BYTE[5][u8(N >> 40)]
      ^ EXTRACT_F16_BYTE[6][u8(N >> 48)] ^ EXTRACT_F16_BYTE[7][u8(N >> 56)];
}
[[gnu::always_inline]] inline u64 embed_f16(u16 n) {
 return EMBED_F16_BYTE[0][u8(n)] ^ EMBED_F16_BYTE[1][u8(n >> 8)];
}

PCLMUL_FN u64 inv_in_f16(u64 N_poly) {
 const u16 n16 = extract_f16(N_poly);
 const u16 inv_n16 = INV16[n16];
 return embed_f16(inv_n16);
}

PCLMUL_FN void inv_batch4(const u64 a[4], u64 out[4]) {
 u64 b1[4], b2[4], b3[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) b1[k] = frob16(a[k]);
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) b2[k] = frob16(b1[k]);
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) b3[k] = frob16(b2[k]);
 u64 beta[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) beta[k] = mul(mul(b1[k], b2[k]), b3[k]);
 u64 N[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) N[k] = mul(a[k], beta[k]);
 u64 N_inv[4];
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) N_inv[k] = inv_in_f16(N[k]);
 #pragma GCC unroll 4
 for (int k = 0; k < 4; ++k) out[k] = mul(beta[k], N_inv[k]);
}

PCLMUL_FN u64 inv_single(u64 a) {
 const u64 b1 = frob16(a);
 const u64 b2 = frob16(b1);
 const u64 b3 = frob16(b2);
 const u64 beta = mul(mul(b1, b2), b3);
 const u64 N = mul(a, beta);
 return mul(beta, inv_in_f16(N));
}

} // namespace gf2_64_pclmul_norm_inv16_compact

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_inv16_compact::inv_batch4;
  using gf2_64_pclmul_norm_inv16_compact::inv_single;
  using gf2_64_pclmul_norm_inv16_compact::init_tables;
  init_tables();
  const size_t T = as.size();
  vector<u64> ans(T);
  size_t i = 0;
  for (; i + 4 <= T; i += 4) {
   u64 b_inv[4];
   inv_batch4(&bs[i], b_inv);
   #pragma GCC unroll 4
   for (int k = 0; k < 4; ++k) ans[i+k] = mul(as[i+k], b_inv[k]);
  }
  for (; i < T; ++i) ans[i] = mul(as[i], inv_single(bs[i]));
  return ans;
 }
};
