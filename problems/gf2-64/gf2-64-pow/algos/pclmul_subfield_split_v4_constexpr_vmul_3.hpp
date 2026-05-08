#pragma once
// pclmul_subfield_split_v4_2 の preprocessing 完全 constexpr 化版:
//
// LN_SIGMA / PW_SIGMA_IDX (各 ~128 KiB) を natural 表現での σ chain で compile-time 構築。
// PCLMUL の constexpr mul は不要 — chain は 16-bit polynomial の (cur<<1) ^ (Q_LOW & -hi)
// で済むので 65535 周しても constexpr step 上限に余裕で収まる。
//
// 利点: init_tables 不要 → 初回 query の cold start 無し、コンパイラの定数畳み込みも狙える。
//
// 実装メモ: subfield_split_constexpr_2.hpp (div 側) と同じ手法。
//   - chain 生成元 BETA を nat 表現の "y" として使い 16-bit shift+XOR でループ。
//   - col[i] = u16(BETA^i) (poly basis) を hardcode、natural→low 変換 byte table を作る。
//   - SUBFIELD_BASIS (subfield 元の low → poly 64bit) は生成元非依存なので v4_2 と同じ。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v4_7 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
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
// natural chain で 65535 周しつつ T_lo/T_hi byte table で nat→low 変換し、両テーブル同時に埋める。
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
VPCLMUL inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec, u64& r0, u64& r1) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i h_idx= _mm256_srli_epi64(prod, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE, indices);
 __m256i result= _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), red_vec);
 r0= _mm256_extract_epi64(result, 0);
 r1= _mm256_extract_epi64(result, 2);
 return result;
}
u64 pow(u64 a, u64 e) {
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

 const u64 N2= mul(a, frob32(a));
 const u16 N= mul(N2, frob16(N2));
 const u64 b= embed_idx(TABLES.PW_SIGMA_IDX[(u32(TABLES.LN_SIGMA[N]) * q) % 65535]);

 // T[i] = a^i for i = 0..15、 binary-tree で 4 層に分けて VPCLMUL 並列化
 u64 T[16]= {1, a, sq(a)};
 // L2: T[3], T[4]
 __m256i T2a= _mm256_set_epi64x(0, T[2], 0, a);
 __m256i T43= mul2(T2a, _mm256_set1_epi64x(T[2]), T[3], T[4]);
 // L3: T[5..8]
 __m256i T4= _mm256_set1_epi64x(T[4]);
 __m256i T56= mul2(T4, T2a, T[5], T[6]);
 mul2(T4, T43, T[7], T[8]);
 // L4: T[9..15] (T[15] は単独 mul)
 __m256i T8= _mm256_set1_epi64x(T[8]);
 mul2(T8, T2a, T[9], T[10]);
 mul2(T8, T43, T[11], T[12]);
 mul2(T8, T56, T[13], T[14]);
 T[15]= mul(T[8], T[7]);
 // メイン loop: 4-bit nibble ごとに frob4 + mul で進める (serial chain なので変更なし)
 int top= 15 - (__builtin_clzll(r) >> 2);
 u64 g= T[(r >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  g= frob4(g);
  u16 chunk= (r >> (4 * i)) & 0xF;
  if(chunk) g= mul(g, T[chunk]);
 }
 return mul(b, g);
}
}  // namespace gf2_64_pow_subfield_split_v4_7
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v4_7::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
