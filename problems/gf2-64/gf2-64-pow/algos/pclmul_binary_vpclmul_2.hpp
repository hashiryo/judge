#pragma once
// pclmul_binary_vpclmul.hpp の小最適化版:
//
// 各 iter で計算する 2 mul はどちらも `a` を掛ける:
//   mul(res, a) と mul(a, a)
// なので b_vec は a の broadcast 1 つで済む。
// _mm256_set_epi64x (4 個の i64 を別々に load) よりも _mm256_set1_epi64x (vpbroadcastq 1 命令)
// のほうが setup が軽い。 1 iter あたり数 cycle のセットアップ削減 (微小)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_pclmul_binary_vpclmul_2 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
// 2 並列 mul (b 共有版): (mul(a0, b), mul(a1, b)) を計算
// b_vec は vpbroadcastq 相当で生成、 a_vec のみ 2 値 set。
VPCLMUL inline void mul2_b_shared(u64 a0, u64 a1, u64 b, u64& r0, u64& r1) {
 __m256i a_vec= _mm256_set_epi64x(0, a1, 0, a0);
 __m256i b_vec= _mm256_set1_epi64x(b);
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i partial= _mm256_xor_si256(prod, red1_shift);
 __m128i p0= _mm256_castsi256_si128(partial);
 __m128i p1= _mm256_extracti128_si256(partial, 1);
 r0= u64(p0[0]) ^ RED[p0[1] >> 60];
 r1= u64(p1[0]) ^ RED[p1[1] >> 60];
}
inline u64 pow_binary(u64 a, u64 e) {
 u64 res= 1;
 while(e) {
  if(e & 1) {
   if(e == 1) {
    res= mul(res, a);
    break;
   }
   // mul(res, a) と mul(a, a) を同じ b=a で 1 vpclmul に
   mul2_b_shared(res, a, a, res, a);
  } else a= sq(a);
  e>>= 1;
 }
 return res;
}
}  // namespace gf2_64_pow_pclmul_binary_vpclmul_2
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pclmul_binary_vpclmul_2::pow_binary;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow_binary(as[i], es[i]);
  return ans;
 }
};
