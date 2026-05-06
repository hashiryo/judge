#pragma once
// PSHUFB (SSSE3) ベースの squaring:
//   GF(2) で a^2 = bit i → bit 2i (cross term ゼロ)。
//   4-bit nibble 単位の spread (nibble 値 v ∈ [0,15] → 8-bit spread, bits at 0,2,4,6) を
//   16 entry テーブルで保持し、PSHUFB で 16 nibble parallel lookup する。
//
// 流れ:
//   1. v = [0:a] を XMM に置く (8 byte が下位)
//   2. lo_n = v & 0x0F   (各 byte の低 nibble)
//   3. hi_n = (v >> 4) & 0x0F (各 byte の高 nibble、SRLI epi16 + AND)
//   4. lo_s = PSHUFB(spread_tbl, lo_n)
//   5. hi_s = PSHUFB(spread_tbl, hi_n)
//   6. squared = unpacklo(lo_s, hi_s) で交互配置 → 128-bit 平方
//   7. squared.hi64 を反 reducible 部分とみなし PDEP 版と同じ shift+XOR で還元
//
// PDEP との比較 (実機):
//   PDEP は Intel ~3 cyc, AMD Zen3 まで ~18 cyc (Zen4 で改善)
//   PSHUFB は arch 共通で 1-1.5 cyc
//   → AMD で有利、Intel/ARM でも互角
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("ssse3")
#endif
#include "../../_shared/_common.hpp"
namespace gf2_64_sq_pshufb {
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
inline u64 sq(u64 a) {
 const __m128i mask_lo= _mm_set1_epi8(0x0F);
 // spread_tbl[v] = bit i ∈ v を bit 2i に展開した 8-bit 値
 //   v=0x0:0x00, 0x1:0x01, 0x2:0x04, 0x3:0x05, 0x4:0x10, ... , 0xF:0x55
 const __m128i spread_tbl= _mm_setr_epi8(
  0x00, 0x01, 0x04, 0x05, 0x10, 0x11, 0x14, 0x15,
  0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55);
 __m128i v= _mm_cvtsi64_si128((i64)a);
 __m128i lo_n= _mm_and_si128(v, mask_lo);
 __m128i hi_n= _mm_and_si128(_mm_srli_epi16(v, 4), mask_lo);
 __m128i lo_s= _mm_shuffle_epi8(spread_tbl, lo_n);
 __m128i hi_s= _mm_shuffle_epi8(spread_tbl, hi_n);
 __m128i squared= _mm_unpacklo_epi8(lo_s, hi_s);
 // squared の下位 64 bit が a の下位 32 bit の平方、上位 64 bit が a の上位 32 bit の平方
 u64 lo= (u64)_mm_cvtsi128_si64(squared);
 u64 h= (u64)_mm_cvtsi128_si64(_mm_srli_si128(squared, 8));
 u64 d= h ^ (h << 1);
 return lo ^ RED[h >> 60] ^ d ^ (d << 3);
}
}  // namespace gf2_64_sq_pshufb
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_pshufb::sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sq(as[i]);
  return ans;
 }
};
