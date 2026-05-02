#pragma once
// Norm-based inv で F_{2^16} subfield を log/exp テーブル直引きする版。
//
// pclmul_norm.hpp は subfield inv に Itoh-Tsujii (6 muls + 15 sqs ≈ 60 cycles) を
// 使うが、F_{2^16} の log/exp テーブル (256 KiB) を持てば 2 lookup (~20 cycles) で済む。
//
// 鍵: nim 基底では F_{2^16} subfield = low 16 bits。既存の basis_change.hpp で
//     poly ↔ nim を 8 byte lookup × 2 で行えるので、subfield inv 全体は
//       poly_to_nim → low 16 bits → log/exp inv → upper 0 → nim_to_poly
//     合計 8 + 2 + 8 = 18 byte lookups (~50 cycles)。
//
// log/exp テーブルは Nimber.hpp と同じ recurrence で構築 (PW[i] = nim 表現での s^i)。
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

namespace gf2_64_pclmul_norm_logexp {
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

// Frobenius (α → α^{2^16}) byte table
inline u64 FROB16_BYTE[8][256];
// F_{2^16} log/exp テーブル (Nimber 互換, nim 基底での値が index/値となる)
inline u16 PW16[65536], LN16[65536];
inline bool inited = false;

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 // Frobenius byte table
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
 // F_{2^16} log/exp テーブル (Nimber.hpp と同じ recurrence)。
 // 第 1 段で poly 基底 (s = 2 を生成元、r(s) = s^16 + ...) で PW[i] = s^i
 PW16[0] = PW16[65535] = 1;
 for (int i = 1; i < 65535; ++i) {
  PW16[i] = u16((PW16[i-1] << 1) ^ (0x1681fu & u16(-(PW16[i-1] >= 0x8000u))));
 }
 // 第 2 段で nim 基底に変換し、変換後の値で LN を引けるようにする
 constexpr u16 f2n[16] = {0x0001u, 0x2827u, 0x392bu, 0x8000u, 0x20fdu, 0x4d1du, 0xde4au, 0x0a17u,
                          0x3464u, 0xe3a9u, 0x6d8du, 0x34bcu, 0xa921u, 0xa173u, 0x0ebcu, 0x0e69u};
 for (int i = 1; i < 65535; ++i) {
  u16 x = PW16[i], y = 0;
  for (; x; x &= x - 1) y ^= f2n[__builtin_ctz(x)];
  PW16[i] = y;
  LN16[y] = u16(i);
 }
 LN16[1] = 0;
}

[[gnu::always_inline]] inline u64 frob16(u64 a) {
 return FROB16_BYTE[0][u8(a)]       ^ FROB16_BYTE[1][u8(a >>  8)]
      ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)]
      ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)]
      ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}

// F_{2^16} subfield 元 N (poly 基底 64-bit) の inv を log/exp で計算 → 64-bit poly に戻す
PCLMUL_FN u64 inv_in_f16_logexp(u64 N_poly) {
 const u64 N_nim = gf2_64_basis::poly_to_nim(N_poly);
 const u16 n16 = u16(N_nim);  // nim 基底では subfield 元の上位 48 bit は 0
 const u16 log_n = LN16[n16];
 const u16 inv_log = u16((65535u - u32(log_n)) % 65535u);
 const u16 inv_n16 = PW16[inv_log];
 return gf2_64_basis::nim_to_poly(u64(inv_n16));
}

PCLMUL_FN u64 inv(u64 a) {
 const u64 b1 = frob16(a);
 const u64 b2 = frob16(b1);
 const u64 b3 = frob16(b2);
 const u64 beta = mul(mul(b1, b2), b3);
 const u64 N = mul(a, beta);
 return mul(beta, inv_in_f16_logexp(N));
}

} // namespace gf2_64_pclmul_norm_logexp

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_logexp::inv;
  using gf2_64_pclmul_norm_logexp::init_tables;
  init_tables();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = mul(as[i], inv(bs[i]));
  return ans;
 }
};
