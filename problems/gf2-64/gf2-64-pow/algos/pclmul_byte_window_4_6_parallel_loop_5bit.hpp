#pragma once
// pclmul_byte_window_4_6_parallel_loop の 5-bit window 版:
//
// 2-lane 分割 (e = e_H·2^32 + e_L, 結合 frob32) はそのままに、window を 4 bit → 5 bit へ。
// 32 bit = 5×6 + 2 なので各 lane は 7 chunk (最上位のみ 2 bit) になり、
// loop が init + 7 反復 → init + 6 反復に減る。結合は既存 FROB32_BYTE / FROB5_BYTE で
// 新テーブル不要 (4_2_parallel_loop_5bit と同じ分割設計)。
//
// T[32] の構築は 4_6 と同形の ×a^2 1 本鎖 (mul2 14 本 + T[31] は単独 mul)。
// 深さ 14 と引き換えに vector 再構築ゼロ・幅最小の 4_6 流を貫く。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window_46_parallel_loop_5bit {
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob5;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
VPCLMUL inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec, u64& r0, u64& r1) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i h_idx= _mm256_srli_epi64(prod, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE, indices);
 __m256i result= _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), red_vec);
 r0= _mm256_castsi256_si128(result)[0];
 r1= _mm256_extracti128_si256(result, 1)[0];
 return result;
}
u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 // T[i] = a^i for i = 0..31、 ×a^2 の 1 本 chain で構築 (4_6 と同形)
 u64 T[32]= {1, a, sq(a)};
 __m256i T2= _mm256_set1_epi64x(T[2]);
 __m256i s= _mm256_set_epi64x(0, T[2], 0, a);
 s= mul2(T2, s, T[3], T[4]);
 s= mul2(T2, s, T[5], T[6]);
 s= mul2(T2, s, T[7], T[8]);
 s= mul2(T2, s, T[9], T[10]);
 s= mul2(T2, s, T[11], T[12]);
 s= mul2(T2, s, T[13], T[14]);
 s= mul2(T2, s, T[15], T[16]);
 s= mul2(T2, s, T[17], T[18]);
 s= mul2(T2, s, T[19], T[20]);
 s= mul2(T2, s, T[21], T[22]);
 s= mul2(T2, s, T[23], T[24]);
 s= mul2(T2, s, T[25], T[26]);
 s= mul2(T2, s, T[27], T[28]);
 s= mul2(T2, s, T[29], T[30]);
 T[31]= mul(T[15], T[16]);
 // メイン loop: 前半/後半 32 bit を 5-bit chunk (下位から 5×6 + 上位 2 bit) の 2-lane lockstep で処理
 const u32 el= u32(e), eh= u32(e >> 32);
 u64 acl= T[el >> 30], ach= T[eh >> 30];
 for(int i= 5; i >= 0; --i) mul2(_mm256_set_epi64x(0, frob5(acl), 0, frob5(ach)), _mm256_set_epi64x(0, T[(el >> (5 * i)) & 31], 0, T[(eh >> (5 * i)) & 31]), ach, acl);
 return mul(frob32(ach), acl);
}
}  // namespace gf2_64_pow_byte_window_46_parallel_loop_5bit
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window_46_parallel_loop_5bit::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
