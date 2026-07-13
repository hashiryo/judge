#pragma once
// メイン loop 16-bit × 4-lane 分割版 (byte_window parallel_loop の 4 分割発展):
//
// e = e_3·2^48 + e_2·2^32 + e_1·2^16 + e_0 と 16 bit ずつ 4 分割すると
//   a^e = frob48(a^{e_3}) · frob32(a^{e_2}) · frob16(a^{e_1}) · a^{e_0}
// 4 lane とも 4 nibble 固定 (init + 3 反復) の同型漸化式なので lockstep で回し、
// mul は mul2 (VPCLMUL) ×2、frob4 は XOR 集約を __m256i で畳む frob4_4lane で 4 lane 同時。
// 直列チェーン: 2-lane 版の init + 7 反復 → init + 3 反復 + 結合 (frob ×3 + mul2 + mul)。
//
// - chunk=0 の conditional skip は無し (T[0]=1 の単位元掛け)
// - e_i = 0 の lane は T[0]=1 が素通しになるだけで正しい (分岐不要)
// - 結合の frob16/32/48 は全て既存 byte table (新テーブル不要)
// - frob4_4lane は byte lookup 32 本を GPR で引きつつ XOR の畳み込みを vector 化。
//   lane の抽出/再構築が毎反復入るのが対価 (loop 4 反復ぶんの短縮に見合うかの実験)
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降, dashboard EPYC 7763 で動作)。
#include <array>
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,vpclmulqdq,avx,avx2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window_parallel_loop {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob48;
using gf2_64_pclmul::FROB4_BYTE;
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
inline array<u64, 4> frob4_4lane(u64 a0, u64 a1, u64 a2, u64 a3) {
 __m256i vec= _mm256_set_epi64x(FROB4_BYTE[0][u8(a3)], FROB4_BYTE[0][u8(a2)], FROB4_BYTE[0][u8(a1)], FROB4_BYTE[0][u8(a0)]);
 for(int i= 1; i < 8; ++i) vec= _mm256_xor_si256(vec, _mm256_set_epi64x(FROB4_BYTE[i][u8(a3 >> (8 * i))], FROB4_BYTE[i][u8(a2 >> (8 * i))], FROB4_BYTE[i][u8(a1 >> (8 * i))], FROB4_BYTE[i][u8(a0 >> (8 * i))]));
 return {u64(_mm256_extract_epi64(vec, 0)), u64(_mm256_extract_epi64(vec, 1)), u64(_mm256_extract_epi64(vec, 2)), u64(_mm256_extract_epi64(vec, 3))};
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
 // メイン loop: 16 bit × 4 lane lockstep (各 lane 4 nibble, init + 3 反復)
 const u16 e0= u16(e), e1= u16(e >> 16), e2= u16(e >> 32), e3= u16(e >> 48);
 u64 A0= T[e0 >> 12], A1= T[e1 >> 12], A2= T[e2 >> 12], A3= T[e3 >> 12];
 for(int i= 2; i >= 0; --i) {
  auto [f0, f1, f2, f3]= frob4_4lane(A0, A1, A2, A3);
  tie(A0, A1)= unpack(mul2(_mm256_set_epi64x(0, f1, 0, f0), _mm256_set_epi64x(0, T[(e1 >> (4 * i)) & 0xF], 0, T[(e0 >> (4 * i)) & 0xF])));
  tie(A2, A3)= unpack(mul2(_mm256_set_epi64x(0, f3, 0, f2), _mm256_set_epi64x(0, T[(e3 >> (4 * i)) & 0xF], 0, T[(e2 >> (4 * i)) & 0xF])));
 }
 // 結合: a^e = frob48(A3)·frob32(A2)·frob16(A1)·A0
 tie(A0, A2)= unpack(mul2(_mm256_set_epi64x(0, frob32(A2), 0, A0), _mm256_set_epi64x(0, frob48(A3), 0, frob16(A1))));
 return mul(A0, A2);
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
