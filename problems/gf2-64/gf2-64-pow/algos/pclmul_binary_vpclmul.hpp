#pragma once
// 素朴 binary exponentiation の VPCLMULQDQ 並列化版:
//
// 各 iteration で必要な 2 mul:
//   - sq(a)     = mul(a, a)         (常時)
//   - mul(res, a)                   (e の最下位 bit が立っている時のみ)
// この 2 つは互いに独立 (どちらも前 iter の res, a を読んで、 次 iter の res, a を作る)。
// VPCLMULQDQ で 1 命令にまとめれば、 mul + sq の 2 issue を 1 issue に短縮。
//
// bit=1 iter で 1 vpclmul、 bit=0 iter で 1 sq (vpclmul は片方無駄になるので使わない)。
// 平均 popcount=32 の e で ~32 issue 削減、 dense_e (popcount~60) でより効果的。
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_pclmul_binary_vpclmul {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
// 2 並列 mul: (mul(a0, b0), mul(a1, b1)) を VPCLMULQDQ + 並列 reduction で計算
VPCLMUL inline void mul2(u64 a0, u64 b0, u64 a1, u64 b1, u64& r0, u64& r1) {
 __m256i a_vec= _mm256_set_epi64x(0, (i64)a1, 0, (i64)a0);
 __m256i b_vec= _mm256_set_epi64x(0, (i64)b1, 0, (i64)b0);
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i partial= _mm256_xor_si256(prod, red1_shift);
 __m128i p0= _mm256_castsi256_si128(partial);
 __m128i p1= _mm256_extracti128_si256(partial, 1);
 r0= (u64)p0[0] ^ RED[(u64)p0[1] >> 60];
 r1= (u64)p1[0] ^ RED[(u64)p1[1] >> 60];
}
inline u64 pow_binary(u64 a, u64 e) {
 u64 res= 1;
 while(e) {
  if(e & 1) {
   if(e == 1) {
    // 最終 iter は res の更新だけで OK (sq(a) 不要)
    res= mul(res, a);
    break;
   }
   // mul(res, a) と sq(a) を vpclmul で 1 命令に
   mul2(res, a, a, a, res, a);
  } else a= sq(a);
  e>>= 1;
 }
 return res;
}
}  // namespace gf2_64_pow_pclmul_binary_vpclmul
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pclmul_binary_vpclmul::pow_binary;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow_binary(as[i], es[i]);
  return ans;
 }
};
