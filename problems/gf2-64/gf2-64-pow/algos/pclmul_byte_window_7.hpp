#pragma once
// pclmul_byte_window_5 の frob4 集約版 (6_2 の (q0,q2) 集約を 2-lane 分割に適用):
//
// 分割・結合は 5 と同一 (e = e_H·2^32 + e_L, 8 nibble × 2 lane, init + 7 反復, frob32)。
// 5 は frob4 を GPR byte table で 2 lane 分引いてから set_epi64x で ymm に載せ直していたが、
// この版は XOR 集約を最初から mul2 の operand 配置 (qword 0, 2) の ymm 上で行い、
// frob の結果ベクタをそのまま mul2 に直結する (q1, q3 は clmul に読まれないので不定で良い)。
// GPR XOR 14 本 + 載せ替え (vmovq ×2 + vinserti128) → vector XOR 7 本 + 2 個詰め set ×8。
//
// 位置づけ (ablation):
//   5 vs 7   = frob 集約そのものの価値 (分割数 2-lane 固定)
//   7 vs 6_2 = 分割数 (init+7 vs init+3) の価値 (集約方式は同型)
// 予想: x64 は直列鎖の長さが支配的なので 6_2 には届かず、5 との差分だけ縮む。
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降, dashboard EPYC 7763 で動作)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window_7 {
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::FROB4_BYTE;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
GNU_TARGET("vpclmulqdq") inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_shift= _mm256_srli_si256(_mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3)), 8);
 __m256i indices= _mm256_srli_si256(_mm256_srli_epi64(prod, 60), 8);
 return _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), _mm256_shuffle_epi8(RED_TABLE, indices));
}
// frob4 ×2 lane: XOR 集約を mul2 の operand 配置 (qword 0, 2) で行い、結果を mul2 に直結
inline __m256i frob4_2lane(u64 a0, u64 a1) {
 __m256i vec= _mm256_set_epi64x(0, FROB4_BYTE[0][u8(a1)], 0, FROB4_BYTE[0][u8(a0)]);
 for(int i= 1; i < 8; ++i) vec= _mm256_xor_si256(vec, _mm256_set_epi64x(0, FROB4_BYTE[i][u8(a1 >> (8 * i))], 0, FROB4_BYTE[i][u8(a0 >> (8 * i))]));
 return vec;
}
inline pair<u64, u64> unpack(const __m256i& vec) { return make_pair(u64(_mm256_extract_epi64(vec, 0)), u64(_mm256_extract_epi64(vec, 2))); }
GNU_TARGET("pclmul,vpclmulqdq") u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 // T[i] = a^i for i = 0..15、 binary-tree で 4 層に分けて VPCLMUL 並列化 (5 と同一)
 u64 T[16]= {1, a, sq(a)};
 // L2: T[3], T[4]
 __m256i T12= _mm256_set_epi64x(0, T[2], 0, a);
 __m256i T34= mul2(T12, _mm256_set1_epi64x(T[2]));
 tie(T[3], T[4])= unpack(T34);
 // L3: T[5..8]
 __m256i T4= _mm256_set1_epi64x(T[4]);
 __m256i T56= mul2(T4, T12);
 tie(T[5], T[6])= unpack(T56);
 tie(T[7], T[8])= unpack(mul2(T4, T34));
 // L4: T[9..15] (T[15] は単独 mul)
 __m256i T8= _mm256_set1_epi64x(T[8]);
 tie(T[9], T[10])= unpack(mul2(T8, T12));
 tie(T[11], T[12])= unpack(mul2(T8, T34));
 tie(T[13], T[14])= unpack(mul2(T8, T56));
 T[15]= mul(T[7], T[8]);
 // メイン loop: 前半/後半 32 bit を 2-lane lockstep で処理 (frob4_2lane を mul2 に直結)
 const u32 el= u32(e), eh= u32(e >> 32);
 u64 acl= T[el >> 28], ach= T[eh >> 28];
 for(int i= 6; i >= 0; --i) tie(ach, acl)= unpack(mul2(frob4_2lane(ach, acl), _mm256_set_epi64x(0, T[(el >> (4 * i)) & 0xF], 0, T[(eh >> (4 * i)) & 0xF])));
 return mul(frob32(ach), acl);
}
}  // namespace gf2_64_pow_byte_window_7
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window_7::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
