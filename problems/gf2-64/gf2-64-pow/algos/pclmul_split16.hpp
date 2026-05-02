#pragma once
// 指数 e を 16-bit ごとに 4 chunk に分けて並列処理する pow バリアント。
//
//   α^e = α^{e_0} · (α^{2^{16}})^{e_1} · (α^{2^{32}})^{e_2} · (α^{2^{48}})^{e_3}
//          where e = e_0 + e_1·2^{16} + e_2·2^{32} + e_3·2^{48}, e_i ∈ [0, 65535]
//
// 4 つの sub-pow は全て独立 (異なる base, 異なる exponent) → OoO で並列実行され
// 単一 pow の critical path が短くなる。各 sub-pow は 16-bit binary exp
// (16 sqs + ~8 muls)、frob16 byte table で base を 3 個追加生成。
//
// 期待: 単一 pow の latency は ~1/3 に短縮、ただし throughput は同等程度。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_pow_split16 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;

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
inline bool inited = false;

PCLMUL_FN void init_frob16_table() {
 if (inited) return;
 inited = true;
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

[[gnu::always_inline]] inline u64 frob16(u64 a) {
 return FROB16_BYTE[0][u8(a)]       ^ FROB16_BYTE[1][u8(a >>  8)]
      ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)]
      ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)]
      ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}

// 16-bit exponent の binary exp
[[gnu::always_inline]] inline u64 pow16(u64 a, u32 e) {
 u64 result = 1;
 while (e) {
  if (e & 1) result = mul(result, a);
  a = sq(a);
  e >>= 1;
 }
 return result;
}

PCLMUL_FN u64 pow(u64 a, u64 e) {
 if (e == 0) return 1;
 // e の 16-bit chunks
 const u32 e0 = u32(e)         & 0xFFFF;
 const u32 e1 = u32(e >> 16)   & 0xFFFF;
 const u32 e2 = u32(e >> 32)   & 0xFFFF;
 const u32 e3 = u32(e >> 48)   & 0xFFFF;
 // Frobenius shifted bases (sequential, 各 frob16 が前の依存)
 const u64 a0 = a;
 const u64 a1 = frob16(a0);
 const u64 a2 = frob16(a1);
 const u64 a3 = frob16(a2);
 // 4 並列 sub-pow (完全に独立、OoO が並列化)
 const u64 r0 = pow16(a0, e0);
 const u64 r1 = pow16(a1, e1);
 const u64 r2 = pow16(a2, e2);
 const u64 r3 = pow16(a3, e3);
 // Tree mul (depth log4 = 2)
 return mul(mul(r0, r1), mul(r2, r3));
}

} // namespace gf2_64_pow_split16

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_split16::pow;
  using gf2_64_pow_split16::init_frob16_table;
  init_frob16_table();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = pow(as[i], es[i]);
  return ans;
 }
};
