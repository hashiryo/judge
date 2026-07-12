#pragma once
// pclmul_byte_window_4_2 のメイン loop 2-lane 分割版:
//
// e = e_H·2^32 + e_L と分けると a^e = frob32(a^{e_H}) · a^{e_L}。
// 両 lane とも 8 nibble 固定 (init + 7 反復) の同型漸化式 (frob4 + mul) なので、
// lockstep で回して mul を mul2 (VPCLMUL) で 2 lane 同時に実行する。
// 直列チェーン init + ~15 反復 → init + 7 反復 + (frob32 + mul) に短縮。
//
// - chunk=0 の conditional skip は廃止して無条件 mul2 (T[0]=1 の単位元掛け。
//   ほぼ毎回 taken だった分岐の mispredict ~1/16 も消える)
// - e_H = 0 (e < 2^32) でも T[0]=1 の lane が素通しになるだけで正しい (分岐不要)
// - 結合の frob32 は既存 FROB32_BYTE をそのまま利用 (新テーブル不要)
// - frob4 は GPR byte table のまま両 lane 並列 (load は独立なので latency はスカラ 1 回分)
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降, dashboard EPYC 7763 で動作)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window_parallel_loop {
using gf2_64_pclmul::frob32;
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
 // T[i] = a^i for i = 0..15、 binary-tree で 4 層に分けて VPCLMUL 並列化 (4_2 と同一)
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
 // メイン loop: 前半/後半 32 bit を 2-lane lockstep で処理
 const u32 el= u32(e), eh= u32(e >> 32);
 u64 acl= T[el >> 28], ach= T[eh >> 28];
 for(int i= 6; i >= 0; --i) {
  acl= frob4(acl), ach= frob4(ach);
  mul2(acl, T[(el >> (4 * i)) & 0xF], ach, T[(eh >> (4 * i)) & 0xF], acl, ach);
 }
 return mul(frob32(ach), acl);
}
}  // namespace gf2_64_pow_byte_window_parallel_loop
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window_parallel_loop::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
