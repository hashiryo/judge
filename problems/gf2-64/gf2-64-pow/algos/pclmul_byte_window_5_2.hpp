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
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_shift= _mm256_srli_si256(_mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3)), 8);
 __m256i indices= _mm256_srli_si256(_mm256_srli_epi64(prod, 60), 8);
 return _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), _mm256_shuffle_epi8(RED_TABLE, indices));
}
inline pair<u64, u64> unpack(const __m256i& vec) { return make_pair(u64(_mm256_extract_epi64(vec, 0)), u64(_mm256_extract_epi64(vec, 2))); }
u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 // T[i] = a^i for i = 0..15、 binary-tree で 4 層に分けて VPCLMUL 並列化 (4_2 と同一)
 u64 T[16]= {1, a, sq(a)};
 // L2: T[3], T[4]
 __m256i T12= _mm256_set_epi64x(0, T[2], 0, a);
 __m256i T34= mul2(T12, _mm256_set1_epi64x(T[2]));
 tie(T[3], T[4])= unpack(T34);
 // L3: T[5..8]
 __m256i T4= _mm256_set1_epi64x(T[4]);
 __m256i T56= mul2(T4, T12);
 tie(T[5], T[6])= unpack(T56);
 __m256i T78= mul2(T4, T34);
 tie(T[7], T[8])= unpack(T78);
 __m256i T910= mul2(T4, T56);
 tie(T[9], T[10])= unpack(T910);
 tie(T[11], T[12])= unpack(mul2(T4, T78));
 tie(T[13], T[14])= unpack(mul2(T4, T910));
 T[15]= mul(T[7], T[8]);
 // メイン loop: 前半/後半 32 bit を 2-lane lockstep で処理
 const u32 el= u32(e), eh= u32(e >> 32);
 u64 acl= T[el >> 28], ach= T[eh >> 28];
 for(int i= 6; i >= 0; --i) {
  tie(ach, acl)= unpack(mul2(_mm256_set_epi64x(0, frob4(acl), 0, frob4(ach)), _mm256_set_epi64x(0, T[(el >> (4 * i)) & 0xF], 0, T[(eh >> (4 * i)) & 0xF])));
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
