#pragma once
// pclmul_byte_window_6 の frob4_4lane 配置換え解消版:
//
// 6 は frob4_4lane が packed ymm (q0..q3) で XOR を畳んだ後 extract ×4 で GPR に下ろし、
// loop 側で set_epi64x ×2 により mul2 の operand 配置 (qword 0, 2) へ載せ直していた。
// _mm256_clmulepi64_epi128(a, b, 0) は各 128bit half の下位 qword しか読まないので、
// XOR 集約を最初から (q0, q2) 配置の ymm 2 本 (lane 0,1 用 / lane 2,3 用) で行えば、
// frob の結果ベクタをそのまま mul2 に直結でき、毎反復の extract ×4 + 再構築 ×2 が消える
// (q1, q3 は不定のままで良い)。vpxor は 7 → 14 本に増えるが、table load は同じ 32 本で
// 2 個詰めは vpinsrq 不要 (vmovq + vinserti128) なので詰めコストはほぼ同等。
// mul2 出力側の unpack は残る — 次の frob の byte 抽出 (アドレス計算) は GPR でしか
// できないため。
//
// 分割・結合は 6 と同一:
//   a^e = frob48(a^{e_3}) · frob32(a^{e_2}) · frob16(a^{e_1}) · a^{e_0}
//   (16 bit × 4 lane, 各 lane 4 nibble, init + 3 反復, 結合は frob ×3 + mul2 + mul)
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
namespace gf2_64_pow_byte_window_6_2 {
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
// frob4 ×4 lane: XOR 集約を最初から mul2 の operand 配置 (qword 0, 2) の 2 本で行う
inline pair<__m256i, __m256i> frob4_4lane2(u64 a0, u64 a1, u64 a2, u64 a3) {
 __m256i vA= _mm256_set_epi64x(0, FROB4_BYTE[0][u8(a1)], 0, FROB4_BYTE[0][u8(a0)]);
 __m256i vB= _mm256_set_epi64x(0, FROB4_BYTE[0][u8(a3)], 0, FROB4_BYTE[0][u8(a2)]);
 for(int i= 1; i < 8; ++i) {
  vA= _mm256_xor_si256(vA, _mm256_set_epi64x(0, FROB4_BYTE[i][u8(a1 >> (8 * i))], 0, FROB4_BYTE[i][u8(a0 >> (8 * i))]));
  vB= _mm256_xor_si256(vB, _mm256_set_epi64x(0, FROB4_BYTE[i][u8(a3 >> (8 * i))], 0, FROB4_BYTE[i][u8(a2 >> (8 * i))]));
 }
 return {vA, vB};
}
inline pair<u64, u64> unpack(const __m256i& vec) { return make_pair(u64(_mm256_extract_epi64(vec, 0)), u64(_mm256_extract_epi64(vec, 2))); }
u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 // T[i] = a^i for i = 0..15、 binary-tree で 4 層に分けて VPCLMUL 並列化 (6 と同一)
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
  auto [vA, vB]= frob4_4lane2(A0, A1, A2, A3);
  tie(A0, A1)= unpack(mul2(vA, _mm256_set_epi64x(0, T[(e1 >> (4 * i)) & 0xF], 0, T[(e0 >> (4 * i)) & 0xF])));
  tie(A2, A3)= unpack(mul2(vB, _mm256_set_epi64x(0, T[(e3 >> (4 * i)) & 0xF], 0, T[(e2 >> (4 * i)) & 0xF])));
 }
 // 結合: a^e = frob48(A3)·frob32(A2)·frob16(A1)·A0
 tie(A0, A2)= unpack(mul2(_mm256_set_epi64x(0, frob32(A2), 0, A0), _mm256_set_epi64x(0, frob48(A3), 0, frob16(A1))));
 return mul(A0, A2);
}
}  // namespace gf2_64_pow_byte_window_6_2
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window_6_2::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
