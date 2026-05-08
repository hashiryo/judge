#pragma once
// pclmul_binary_vpclmul_3 + 最終 RED lookup も SIMD 化:
//
// step() の reduction パートは partial を作った後、 各 lane で
//   r_k = lo_k ^ red1_k ^ RED[h_k >> 60]
// と RED 補正だけスカラー残ってた。 これを SIMD で:
//   1. h_idx = srli_epi64(partial, 60)        各 64-bit element を >> 60
//   2. indices = srli_si256(h_idx, 8)         lane 内で byte 8→0 にシフト → byte 0 = h_k>>60
//   3. red_vec = shuffle_epi8(red_table, indices)   PSHUFB で RED lookup
//   4. final = xor(partial, red_vec)          partial に低 byte だけ RED 補正を XOR
//   5. lane k の low 64 を取り出して res, a に代入
//
// red_table は 16-byte/lane で RED[0..7] + 0 padding を 2 lane に複製。
// PSHUFB は 1 lane 内 byte shuffle なので index byte 0 だけ意味あり、 他は 0 (RED[0]=0 で安全)。
//
// 期待: スカラー 2 chain (~5-7 cyc) → SIMD 4 op (~4 cyc) で 1-3 cyc/iter 短縮。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_pclmul_binary_vpclmul_4 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
const __m256i red_table= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
VPCLMUL inline void step(u64& res, u64& a) {
 __m256i a_vec= _mm256_set_epi64x(0, a, 0, res);
 __m256i b_vec= _mm256_set1_epi64x(a);
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i partial= _mm256_xor_si256(prod, red1_shift);
 // RED 補正: PSHUFB で h_k >> 60 を index に lookup → 各 lane 低 byte に RED[h_k>>60]
 __m256i h_idx= _mm256_srli_epi64(partial, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(red_table, indices);
 __m256i final_vec= _mm256_xor_si256(partial, red_vec);
 __m128i p0= _mm256_castsi256_si128(final_vec);
 __m128i p1= _mm256_extracti128_si256(final_vec, 1);
 res= (u64)p0[0];
 a= (u64)p1[0];
}
inline u64 pow_binary(u64 a, u64 e) {
 u64 res= 1;
 while(e) {
  if(e & 1) {
   if(e == 1) {
    res= mul(res, a);
    break;
   }
   step(res, a);
  } else a= sq(a);
  e>>= 1;
 }
 return res;
}
}  // namespace gf2_64_pow_pclmul_binary_vpclmul_4
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pclmul_binary_vpclmul_4::pow_binary;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow_binary(as[i], es[i]);
  return ans;
 }
};
