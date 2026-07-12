#pragma once
// pclmul_byte_window_4_2_parallel_loop の 5-bit window 版:
//
// 2-lane 分割 (e = e_H·2^32 + e_L, 結合 frob32) はそのままに、window を 4 bit → 5 bit へ。
// 32 bit = 5×6 + 2 なので各 lane は 7 chunk (最上位のみ 2 bit) になり、
// loop が init + 7 反復 → init + 6 反復に減る。分割点が 32 bit のままなので
// 結合は既存 FROB32_BYTE、chunk 側も既存 FROB5_BYTE で新テーブル不要。
//
// トレードオフ (実測枠):
//   + 直列 loop が 1 反復 (frob + mul2) 減る
//   - T が a^0..a^15 → a^0..a^31 になり、構築 tree が 6 mul2 → 14 mul2 (+1 層) に増える
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降, dashboard EPYC 7763 で動作)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window_parallel_loop_5bit {
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob5;
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
 // T[i] = a^i for i = 0..31、 binary-tree で 5 層に分けて VPCLMUL 並列化
 u64 T[32]= {1, a, sq(a)};
 // L2: T[3], T[4]
 mul2(T[2], a, T[2], T[2], T[3], T[4]);
 // L3: T[5..8]
 mul2(T[4], a, T[4], T[2], T[5], T[6]);
 mul2(T[3], T[4], T[4], T[4], T[7], T[8]);
 // L4: T[9..16]
 mul2(T[8], a, T[8], T[2], T[9], T[10]);
 mul2(T[8], T[3], T[8], T[4], T[11], T[12]);
 mul2(T[8], T[5], T[8], T[6], T[13], T[14]);
 mul2(T[8], T[7], T[8], T[8], T[15], T[16]);
 // L5: T[17..31] (T[31] は単独 mul)
 mul2(T[16], a, T[16], T[2], T[17], T[18]);
 mul2(T[16], T[3], T[16], T[4], T[19], T[20]);
 mul2(T[16], T[5], T[16], T[6], T[21], T[22]);
 mul2(T[16], T[7], T[16], T[8], T[23], T[24]);
 mul2(T[16], T[9], T[16], T[10], T[25], T[26]);
 mul2(T[16], T[11], T[16], T[12], T[27], T[28]);
 mul2(T[16], T[13], T[16], T[14], T[29], T[30]);
 T[31]= mul(T[16], T[15]);
 // メイン loop: 前半/後半 32 bit を 5-bit chunk (下位から 5×6 + 上位 2 bit) の 2-lane lockstep で処理
 const u32 el= u32(e), eh= u32(e >> 32);
 u64 acl= T[el >> 30], ach= T[eh >> 30];
 for(int i= 5; i >= 0; --i) {
  acl= frob5(acl), ach= frob5(ach);
  mul2(acl, T[(el >> (5 * i)) & 31], ach, T[(eh >> (5 * i)) & 31], acl, ach);
 }
 return mul(frob32(ach), acl);
}
}  // namespace gf2_64_pow_byte_window_parallel_loop_5bit
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window_parallel_loop_5bit::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
