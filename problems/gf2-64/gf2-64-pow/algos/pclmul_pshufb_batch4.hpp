#pragma once
// PSHUFB-based SIMD batch4 pow:
//
//   4 個の (a, e) ペアをまとめて処理する。各 iteration で:
//     - SIMD で 4 個の squaring を PSHUFB nibble-lookup で並列実行
//     - 4 個の mul は scalar PCLMUL で順次 (PCLMUL に SIMD 版がない)
//
// PSHUFB squaring の仕組み:
//   GF(2) で y → y^2 は "bit i → bit 2i" の spread (cross term ゼロ)。
//   入力 byte b の低 nibble (4 bit) と高 nibble (4 bit) を、それぞれ
//   16 エントリ table で 8-bit 出力に展開:
//     spread_table[i] = i (4-bit) を 1-bit gap で spread した 8-bit
//   2 PSHUFB + unpack で 64-bit → 128-bit spread 完了。
//
//   AVX2 _mm_shuffle_epi8 は 128-bit (16 byte) 並列 lookup。
//   2 個の u64 を __m128i に詰めれば 1 PSHUFB 呼びで 2 値の low-nibble lookup
//   が並列に完了 → 2-way SIMD squaring。
//
// Note: scalar PDEP-sq でも 1 cycle。PSHUFB sq は throughput 改善期待だが、
//       Zen 3 で PDEP 速いので gain は小さいかもしれない。ARM NEON では
//       PSHUFB 相当 (vqtbl1q_u8) があるので意味がより大きい (PDEP 無いから)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2,sse4.1,ssse3")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2,sse4.1,ssse3")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_pow_pshufb_batch4 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;

// 4-bit nibble → 8-bit spread (bit i → bit 2i within byte)
//   0b0000 → 0b00000000 = 0x00
//   0b0001 → 0b00000001 = 0x01
//   0b0010 → 0b00000100 = 0x04
//   0b0011 → 0b00000101 = 0x05
//   0b1111 → 0b01010101 = 0x55
inline constexpr u8 SPREAD_TABLE[16] = {
 0x00, 0x01, 0x04, 0x05, 0x10, 0x11, 0x14, 0x15,
 0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55,
};

// 2 個の u64 を SIMD で squaring (spread 部分のみ、reduce は scalar)
PCLMUL_FN void sq_pair_simd(u64& v0, u64& v1) {
 __m128i tbl = _mm_loadu_si128(reinterpret_cast<const __m128i*>(SPREAD_TABLE));
 __m128i lo_mask = _mm_set1_epi8(0x0F);
 __m128i input = _mm_set_epi64x((long long) v1, (long long) v0);
 __m128i lo_nibs = _mm_and_si128(input, lo_mask);
 __m128i hi_nibs = _mm_and_si128(_mm_srli_epi64(input, 4), lo_mask);
 __m128i sp_lo = _mm_shuffle_epi8(tbl, lo_nibs);
 __m128i sp_hi = _mm_shuffle_epi8(tbl, hi_nibs);
 // unpacklo: byte 2i = sp_lo[i], byte 2i+1 = sp_hi[i] for i = 0..7 → v0 spread (128-bit)
 // unpackhi: 同 for i = 8..15 → v1 spread (128-bit)
 __m128i spread_v0 = _mm_unpacklo_epi8(sp_lo, sp_hi);
 __m128i spread_v1 = _mm_unpackhi_epi8(sp_lo, sp_hi);
 v0 = reduce(spread_v0);
 v1 = reduce(spread_v1);
}

// 4 個まとめて SIMD squaring
PCLMUL_FN void sq_quad_simd(u64& v0, u64& v1, u64& v2, u64& v3) {
 sq_pair_simd(v0, v1);
 sq_pair_simd(v2, v3);
}

// 4-way batched binary exponentiation
PCLMUL_FN void pow_batch4(const u64 a[4], const u64 e[4], u64 out[4]) {
 u64 acc0 = 1, acc1 = 1, acc2 = 1, acc3 = 1;
 u64 b0 = a[0], b1 = a[1], b2 = a[2], b3 = a[3];
 const u64 e0 = e[0], e1 = e[1], e2 = e[2], e3 = e[3];
 const u64 e_or = e0 | e1 | e2 | e3;
 if (e_or == 0) { out[0] = out[1] = out[2] = out[3] = 1; return; }
 const int top_bit = 63 - __builtin_clzll(e_or);
 for (int bit = 0; bit <= top_bit; ++bit) {
  if ((e0 >> bit) & 1) acc0 = mul(acc0, b0);
  if ((e1 >> bit) & 1) acc1 = mul(acc1, b1);
  if ((e2 >> bit) & 1) acc2 = mul(acc2, b2);
  if ((e3 >> bit) & 1) acc3 = mul(acc3, b3);
  if (bit < top_bit) sq_quad_simd(b0, b1, b2, b3);
 }
 out[0] = acc0; out[1] = acc1; out[2] = acc2; out[3] = acc3;
}

PCLMUL_FN u64 pow_single(u64 a, u64 e) {
 // Fall-back for tail when T not divisible by 4
 if (e == 0) return 1;
 u64 result = 1;
 while (e) {
  if (e & 1) result = mul(result, a);
  // PSHUFB sq, 1-way (use sq_pair_simd with dummy)
  u64 dummy = 0;
  sq_pair_simd(a, dummy);
  e >>= 1;
 }
 return result;
}

} // namespace gf2_64_pow_pshufb_batch4

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pshufb_batch4::pow_batch4;
  using gf2_64_pow_pshufb_batch4::pow_single;
  const size_t T = as.size();
  vector<u64> ans(T);
  size_t i = 0;
  for (; i + 4 <= T; i += 4) {
   u64 a4[4] = {as[i], as[i+1], as[i+2], as[i+3]};
   u64 e4[4] = {es[i], es[i+1], es[i+2], es[i+3]};
   u64 r4[4];
   pow_batch4(a4, e4, r4);
   ans[i] = r4[0]; ans[i+1] = r4[1]; ans[i+2] = r4[2]; ans[i+3] = r4[3];
  }
  for (; i < T; ++i) ans[i] = pow_single(as[i], es[i]);
  return ans;
 }
};
