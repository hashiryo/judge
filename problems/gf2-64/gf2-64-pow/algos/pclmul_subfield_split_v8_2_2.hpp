#pragma once
// pclmul_subfield_split_v8 の frob4_4lane 配置換え解消版 (byte_window_6_2 の subfield 版):
//
// v8 は frob4_4lane が packed ymm (q0..q3) で XOR を畳んだ後 extract ×4 で GPR に下ろし、
// loop 側で set_epi64x ×2 により mul2 の operand 配置 (qword 0, 2) へ載せ直していた。
// XOR 集約を最初から (q0, q2) 配置の ymm 2 本 (lane 0,1 用 / lane 2,3 用) で行い、
// frob の結果ベクタをそのまま mul2 に直結して毎反復の extract ×4 + 再構築 ×2 を消す
// (q1, q3 は clmul に読まれないので不定で良い)。vpxor は 7 → 14 本に増えるが、
// table load は同じ 32 本で 2 個詰めは vpinsrq 不要なので詰めコストはほぼ同等。
//
// 分割・結合は v8 と同一:
//   a^r = frob36(a^{r_3}) · frob24(a^{r_2}) · frob12(a^{r_1}) · a^{r_0}
//   (12 bit × 4 lane, 各 lane 3 nibble, init + 2 反復, b は lane0 最下位に織り込み)
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v8_2 {
using gf2_64_pclmul::frob12;
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob24;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob36;
using gf2_64_pclmul::FROB4_BYTE;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
// embed: low 16-bit subfield 識別子 → 64-bit poly 表現 (subfield 元の埋め込み)
constexpr u64 embed_idx(u16 idx) {
 static constexpr auto EMBED= []() {
  u64 SUBFIELD_BASIS[]= {0x0000000000000001ULL, 0x5fbfaec6aeac0002ULL, 0xb06c601895640004ULL, 0xb013b5277b7c0008ULL, 0xb5ebb915248a0010ULL, 0x109bb25b2c600020ULL, 0xbf3bd95bd4190040ULL, 0x0fc66342279b0080ULL, 0xb6418f5e57c50100ULL, 0xaa194bd4b83f0200ULL, 0x1b5217b4dcc70400ULL, 0xbb06fa73867a0800ULL, 0x006fd55b23331000ULL, 0x4ae8fb39198c2000ULL, 0xfbd141b29b4f4000ULL, 0x1d9ce1776be78000ULL};
  array<array<u64, 256>, 2> t{};
  for(int half= 0; half < 2; ++half)
   for(int i= 0; i < 256; ++i) {
    u64 v= 0;
    for(int b= 0; b < 8; ++b)
     if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
    t[half][i]= v;
   }
  return t;
 }();
 return EMBED[0][u8(idx)] ^ EMBED[1][idx >> 8];
}
// LN_SIGMA[u16(BETA^k)] = k, PW_SIGMA_IDX[k] = u16(BETA^k)
struct Tables {
 u16 LN_SIGMA[65536];
 u16 PW_SIGMA_IDX[65535];
};
constexpr auto TABLES= []() {
 // col[i] = u16(BETA^i) — BETA = 0x1f1af3ec55a22e02 の poly basis 累乗の低 16 bit
 u16 col[]= {1U, 11778U, 7028U, 51115U, 48663U, 26081U, 17458U, 40223U, 30334U, 42368U, 14380U, 2223U, 49688U, 11217U, 44239U, 63445U};
 u16 T_lo[256]= {}, T_hi[256]= {};
 for(int v= 0; v < 256; ++v) {
  u16 lo= 0, hi= 0;
  for(int j= 0; j < 8; ++j)
   if((v >> j) & 1) {
    lo^= col[j];
    hi^= col[j + 8];
   }
  T_lo[v]= lo;
  T_hi[v]= hi;
 }
 Tables t{};
 // BETA の最小多項式は y^16 + y^5 + y^3 + y^2 + 1, lower bits = 0x002D
 u16 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u16 lo= T_lo[u8(cur)] ^ T_hi[cur >> 8];
  t.LN_SIGMA[lo]= u16(k);
  t.PW_SIGMA_IDX[k]= lo;
  cur= u16(cur << 1) ^ (0x002DU & -u16(cur >> 15));
 }
 t.LN_SIGMA[0]= 0;
 return t;
}();
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
GNU_TARGET("vpclmulqdq") inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec) {
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
GNU_TARGET("pclmul,vpclmulqdq") u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 if(!a) return 0;
 constexpr u64 M_VAL= (~u64(0)) / 65535u;
 const u16 q= e / M_VAL;
 const u64 r= e - M_VAL * q;
 if(!r) {
  const u64 N2= mul(a, frob32(a));
  const u16 N= mul(N2, frob16(N2));
  return embed_idx(TABLES.PW_SIGMA_IDX[(u32(TABLES.LN_SIGMA[N]) * q) % 65535]);
 }

 // T[i] = a^i for i = 0..15、 binary-tree で 4 層に分けて VPCLMUL 並列化
 u64 T[16]= {1, a};
 u64 B;
 tie(T[2], B)= unpack(mul2(_mm256_set1_epi64x(a), _mm256_set_epi64x(0, frob32(a), 0, a)));

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
 tie(T[15], B)= unpack(mul2(_mm256_set_epi64x(0, frob16(B), 0, T[7]), _mm256_set_epi64x(0, B, 0, T[8])));
 // メイン loop 4-lane 化: 各 lane 3 nibble, init + 2 反復 (r ≥ 2^48 のときだけ +1)
 const u32 r0= r & 0xFFF, r1= (r >> 12) & 0xFFF, r2= (r >> 24) & 0xFFF, r3= u32(r >> 36);
 const u64 bt= mul(embed_idx(TABLES.PW_SIGMA_IDX[(u32(TABLES.LN_SIGMA[u16(B)]) * q) % 65535]), T[r0 & 0xF]);
 const int it= int(r >> 48) + 2;
 u64 A0= T[(r0 >> (4 * it)) & 0xF], A1= T[(r1 >> (4 * it)) & 0xF], A2= T[(r2 >> (4 * it)) & 0xF], A3= T[(r3 >> (4 * it)) & 0xF];
 for(int i= it - 1; i >= 1; --i) {
  auto [vA, vB]= frob4_4lane2(A0, A1, A2, A3);
  tie(A0, A1)= unpack(mul2(vA, _mm256_set_epi64x(0, T[(r1 >> (4 * i)) & 0xF], 0, T[(r0 >> (4 * i)) & 0xF])));
  tie(A2, A3)= unpack(mul2(vB, _mm256_set_epi64x(0, T[(r3 >> (4 * i)) & 0xF], 0, T[(r2 >> (4 * i)) & 0xF])));
 }
 // 最終段 (i=0): lane0 に bt (= b·T[r&15]) を織り込み
 auto [vA, vB]= frob4_4lane2(A0, A1, A2, A3);
 tie(A0, A1)= unpack(mul2(vA, _mm256_set_epi64x(0, T[r1 & 0xF], 0, bt)));
 tie(A2, A3)= unpack(mul2(vB, _mm256_set_epi64x(0, T[r3 & 0xF], 0, T[r2 & 0xF])));
 // 結合: a^e = frob36(A3)·frob24(A2)·frob12(A1)·A0 (A0 は b 織り込み済み)
 tie(A0, A2)= unpack(mul2(_mm256_set_epi64x(0, frob24(A2), 0, A0), _mm256_set_epi64x(0, frob36(A3), 0, frob12(A1))));
 return mul(A0, A2);
}
}  // namespace gf2_64_pow_subfield_split_v8_2
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v8_2::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
