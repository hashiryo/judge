#pragma once
// pclmul_binary_vpclmul_2.hpp の更なる小整理:
//
// binary pow の各 set bit iter で必要な計算は常に:
//   res ← mul(res, a)
//   a   ← mul(a, a)
// 入力は (res, a) のみ、 出力も同じ 2 変数 → 関数は 2 引数の参照で表せる。
// b 共有 broadcast はそのまま (b = a)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_pclmul_binary_vpclmul_3 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
// res ← mul(res, a), a ← mul(a, a) を vpclmul で同時に実行
VPCLMUL inline void step(u64& res, u64& a) {
 __m256i a_vec= _mm256_set_epi64x(0, a, 0, res);  // [_, a, _, res]
 __m256i b_vec= _mm256_set1_epi64x(a);            // a を全 lane broadcast
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i partial= _mm256_xor_si256(prod, red1_shift);
 __m128i p0= _mm256_castsi256_si128(partial);
 __m128i p1= _mm256_extracti128_si256(partial, 1);
 res= u64(p0[0]) ^ RED[p0[1] >> 60];
 a= u64(p1[0]) ^ RED[p1[1] >> 60];
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
}  // namespace gf2_64_pow_pclmul_binary_vpclmul_3
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pclmul_binary_vpclmul_3::pow_binary;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow_binary(as[i], es[i]);
  return ans;
 }
};
