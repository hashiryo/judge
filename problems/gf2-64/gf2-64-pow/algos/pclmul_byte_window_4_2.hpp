#pragma once
// pclmul_byte_window_3 + VPCLMULQDQ で T[2..15] precompute を並列化:
//
// 元の sequential T[i] = T[i-1]·a (14 直列 mul) を binary tree で組み直し、
// 独立な mul を 2 個ずつ VPCLMUL で同時実行する。 critical path 14 → 4 layer まで短縮。
//
// 構造:
//   L0: T[1] = a
//   L1: T[2] = a·a                                    (1 mul)
//   L2: T[3] = T[2]·a, T[4] = T[2]·T[2]              (1 vpclmul)
//   L3: T[5] = T[4]·a, T[6] = T[4]·T[2]              (1 vpclmul)
//       T[7] = T[3]·T[4], T[8] = T[4]·T[4]           (1 vpclmul)
//   L4: T[9] = T[8]·a,  T[10] = T[8]·T[2]            (1 vpclmul)
//       T[11] = T[8]·T[3], T[12] = T[8]·T[4]         (1 vpclmul)
//       T[13] = T[8]·T[5], T[14] = T[8]·T[6]         (1 vpclmul)
//       T[15] = T[8]·T[7]                            (1 mul)
//
// メイン loop (acc に対する frob4 + conditional mul) は serial dependency なので変更なし。
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降, dashboard EPYC 7763 で動作)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window_vpclmul {
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
// 2 並列 mul: (mul(a0, b0), mul(a1, b1)) を VPCLMULQDQ + 並列 reduction で計算
VPCLMUL inline void mul2(u64 a0, u64 b0, u64 a1, u64 b1, u64& r0, u64& r1) {
 __m256i a_vec= _mm256_set_epi64x(0, a1, 0, a0);
 __m256i b_vec= _mm256_set_epi64x(0, b1, 0, b0);
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
u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 // T[i] = a^i for i = 0..15、 binary-tree で 4 層に分けて VPCLMUL 並列化
 u64 T[16]= {1, a, sq(a)};
 // L2: T[3], T[4]
 mul2(T[2], a, T[2], T[2], T[3], T[4]);
 // L3: T[5..8]
 mul2(T[4], a, T[4], T[2], T[5], T[6]);
 mul2(T[3], T[4], T[4], T[4], T[7], T[8]);
 // L4: T[9..15] (T[15] は単独 mul)
 mul2(T[8], a, T[8], T[2], T[9], T[10]);
 mul2(T[8], T[3], T[8], T[4], T[11], T[12]);
 mul2(T[8], T[5], T[8], T[6], T[13], T[14]);
 T[15]= mul(T[8], T[7]);
 // メイン loop: 4-bit nibble ごとに frob4 + mul で進める (serial chain なので変更なし)
 int top= 15 - (__builtin_clzll(e) >> 2);
 u64 acc= T[(e >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob4(acc);
  u16 chunk= (e >> (4 * i)) & 0xF;
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
}  // namespace gf2_64_pow_byte_window_vpclmul
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window_vpclmul::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
