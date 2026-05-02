#pragma once
// PCLMUL + PDEP 二乗 + 4-bit sliding window exponentiation。
//
// 標準 binary exp: ~32 muls + 64 sqs。
// 4-bit window: precompute {a, a^3, ..., a^{15}} (= 7 muls)、e の 1 のクラスタ
// (最大 4 bit 連続) を 1 mul で消化 → 平均 ~16 muls + 64 sqs。
// 大局的には mul 数を半減させる狙い。
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

namespace gf2_64_pow_pdep_window {
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

PCLMUL_FN u64 pow(u64 a, u64 e) {
 if (e == 0) return 1;
 // odd 冪のみ precompute: tab[k] = a^{2k+1} for k = 0..7 → a, a^3, ..., a^{15}
 u64 tab[8];
 tab[0] = a;
 const u64 a2 = sq(a);
 for (int k = 1; k < 8; ++k) tab[k] = mul(tab[k-1], a2);

 int msb = 63 - __builtin_clzll(e);
 u64 res = 1;
 int i = msb;
 while (i >= 0) {
  if (!((e >> i) & 1)) {
   res = sq(res);
   --i;
   continue;
  }
  // 最大 4 bit window を取り出し、奇数化
  int w = (i >= 3) ? 4 : (i + 1);
  unsigned win = unsigned((e >> (i - w + 1)) & ((1u << w) - 1));
  int trail = __builtin_ctz(win);
  win >>= trail;
  int actual_w = w - trail;
  for (int s = 0; s < actual_w; ++s) res = sq(res);
  res = mul(res, tab[(win - 1) >> 1]);
  i -= w;
  for (int s = 0; s < trail; ++s) res = sq(res);
 }
 return res;
}

} // namespace gf2_64_pow_pdep_window

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pdep_window::pow;
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = pow(as[i], es[i]);
  return ans;
 }
};
