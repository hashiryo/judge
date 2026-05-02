#pragma once
// Frobenius basis precompute pow:
//
//   α^e = ⊕(× over set bits k of e)  α^{2^k}
//
// 各 α^{2^k} (k = 0..63) を per-call で precompute すれば、e の bit ごとに
// その値を tree mul するだけで pow 完了。
//
// 工夫: α^{2^k} の chain は 4 stream に分割して並列計算
//   stream 0: α, sq(α), sq(sq(α)), ..., α^{2^{15}}      (16 sqs)
//   stream 1: frob16(α), sq(...), ..., α^{2^{31}}        (16 sqs)
//   stream 2: frob16(frob16(α)), sq(...), ..., α^{2^{47}} (16 sqs)
//   stream 3: frob16(frob16(frob16(α))), sq(...), ..., α^{2^{63}} (16 sqs)
// → 4 stream の sq が独立、OoO で並列実行 (実効 ~16 sqs latency)
//
// pclmul_split16 との違い:
//   split16 は 4 つの「サブ pow」を独立に走らせる
//   こちらは α^{2^k} を全部 precompute してから tree mul
//   → 後者の方が tree mul の depth が log2(popcount(e)) と短い (= 5 levels)
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

namespace gf2_64_pow_frob_basis {
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

PCLMUL_FN u64 pow(u64 a, u64 e) {
 if (e == 0) return 1;
 // Precompute α^{2^k} for k = 0..63 を 4 stream 並列で生成
 u64 frob[64];
 frob[0]  = a;
 frob[16] = frob16(a);
 frob[32] = frob16(frob[16]);
 frob[48] = frob16(frob[32]);
 // 4 stream の sq を interleave (compiler/OoO が並列発行する)
 #pragma GCC unroll 15
 for (int i = 1; i < 16; ++i) {
  frob[i]      = sq(frob[i - 1]);
  frob[16 + i] = sq(frob[16 + i - 1]);
  frob[32 + i] = sq(frob[32 + i - 1]);
  frob[48 + i] = sq(frob[48 + i - 1]);
 }
 // Tree mul over selected bits of e
 // popcount(e) 個を 2 つずつ pairwise 掛ける → 深さ log2(popcount) の木
 u64 selected[64];
 int n = 0;
 for (int k = 0; k < 64; ++k) {
  if ((e >> k) & 1) selected[n++] = frob[k];
 }
 // Tree reduction (in-place)
 while (n > 1) {
  int new_n = 0;
  int i = 0;
  for (; i + 1 < n; i += 2) {
   selected[new_n++] = mul(selected[i], selected[i + 1]);
  }
  if (i < n) selected[new_n++] = selected[i];
  n = new_n;
 }
 return selected[0];
}

} // namespace gf2_64_pow_frob_basis

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_frob_basis::pow;
  using gf2_64_pow_frob_basis::init_frob16_table;
  init_frob16_table();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = pow(as[i], es[i]);
  return ans;
 }
};
