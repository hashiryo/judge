#pragma once
// pclmul_byte_window_6_2 の「両 qword lane 活用」版:
//
// _mm256_clmulepi64_epi128 の imm8 は bit0 で a 側、bit4 で b 側の qword (low/high) を
// 選ぶので、imm=0x00 (low×low) と imm=0x11 (high×high) の 2 種類の mul2 を用意すると、
// 1 本の ymm の q0..q3 全部に lane を詰めたまま 4 lane 分の乗算ができる:
//   facc = (f0, f1, f2, f3), tvec = (T0, T1, T2, T3) に対し
//   mul2  (imm 0x00): half0 = f0×T0, half1 = f2×T2 → (A0, A2)
//   mul2h (imm 0x11): half0 = f1×T1, half1 = f3×T3 → (A1, A3)
// reduction は prod の配置が imm に依らず同じなので完全に共通 (template で imm だけ変更)。
//
// 6_2 との違い:
//   - frob4 の XOR 集約が packed 1 本 (vpxor 7 本) に戻り、かつ 6 にあった
//     extract ×4 も不要 (mul2/mul2h が q0..q3 を直接消費する)
//   - T operand の詰めも 2 本 → packed 1 本
//   - 遊んでいた q1, q3 lane が働くので live な ymm が半分に
//
// 分割・結合は 6 系と同一:
//   a^e = frob48(a^{e_3}) · frob32(a^{e_2}) · frob16(a^{e_1}) · a^{e_0}
//   (16 bit × 4 lane, 各 lane 4 nibble, init + 3 反復, 結合は frob ×3 + mul2 + mul)
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降, dashboard EPYC 7763 で動作)。
#include <array>
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window_6_3 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob48;
using gf2_64_pclmul::FROB4_BYTE;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
// 2 並列 mul: IMM で各 128bit half のどの qword を掛けるか選ぶ (0x00: low×low, 0x11: high×high)。
// reduction は prod 配置が IMM に依らないので共通。
template <int IMM> GNU_TARGET("vpclmulqdq") inline __m256i mul2t(const __m256i& a_vec, const __m256i& b_vec) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, IMM);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_shift= _mm256_srli_si256(_mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3)), 8);
 __m256i indices= _mm256_srli_si256(_mm256_srli_epi64(prod, 60), 8);
 return _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), _mm256_shuffle_epi8(RED_TABLE, indices));
}
GNU_TARGET("vpclmulqdq") inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec) { return mul2t<0x00>(a_vec, b_vec); }
GNU_TARGET("vpclmulqdq") inline __m256i mul2h(const __m256i& a_vec, const __m256i& b_vec) { return mul2t<0x11>(a_vec, b_vec); }
// frob4 ×4 lane: packed (q0..q3) のまま XOR 集約し、vector のまま返す (extract 不要)
inline __m256i frob4_4lane(u64 a0, u64 a1, u64 a2, u64 a3) {
 __m256i vec= _mm256_set_epi64x(FROB4_BYTE[0][u8(a3)], FROB4_BYTE[0][u8(a2)], FROB4_BYTE[0][u8(a1)], FROB4_BYTE[0][u8(a0)]);
 for(int i= 1; i < 8; ++i) vec= _mm256_xor_si256(vec, _mm256_set_epi64x(FROB4_BYTE[i][u8(a3 >> (8 * i))], FROB4_BYTE[i][u8(a2 >> (8 * i))], FROB4_BYTE[i][u8(a1 >> (8 * i))], FROB4_BYTE[i][u8(a0 >> (8 * i))]));
 return vec;
}
inline pair<u64, u64> unpack(const __m256i& vec) { return make_pair(u64(_mm256_extract_epi64(vec, 0)), u64(_mm256_extract_epi64(vec, 2))); }
GNU_TARGET("pclmul,vpclmulqdq") u64 pow(u64 a, u64 e) {
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
 // メイン loop: 16 bit × 4 lane lockstep (各 lane 4 nibble, init + 3 反復)。
 // facc/tvec は packed のまま、mul2 (low×low) と mul2h (high×high) で 4 lane 同時に掛ける。
 const u16 e0= u16(e), e1= u16(e >> 16), e2= u16(e >> 32), e3= u16(e >> 48);
 u64 A0= T[e0 >> 12], A1= T[e1 >> 12], A2= T[e2 >> 12], A3= T[e3 >> 12];
 for(int i= 2; i >= 0; --i) {
  __m256i facc= frob4_4lane(A0, A1, A2, A3);
  __m256i tvec= _mm256_set_epi64x(T[(e3 >> (4 * i)) & 0xF], T[(e2 >> (4 * i)) & 0xF], T[(e1 >> (4 * i)) & 0xF], T[(e0 >> (4 * i)) & 0xF]);
  tie(A0, A2)= unpack(mul2(facc, tvec));
  tie(A1, A3)= unpack(mul2h(facc, tvec));
 }
 // 結合: a^e = frob48(A3)·frob32(A2)·frob16(A1)·A0
 tie(A0, A2)= unpack(mul2(_mm256_set_epi64x(0, frob32(A2), 0, A0), _mm256_set_epi64x(0, frob48(A3), 0, frob16(A1))));
 return mul(A0, A2);
}
}  // namespace gf2_64_pow_byte_window_6_3
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window_6_3::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
