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
//
// main loop 2-lane 分割 (この版の追加分):
//   r = r_H·2^24 + r_L と分けて a^r = frob24(a^{r_H})·a^{r_L}。両 lane を lockstep で回して
//   mul を mul2s (VPCLMUL, GPR in/out) で 2 lane 同時実行。frob4 は GPR byte table のまま
//   両 lane 並列 (load 独立なので latency はスカラ 1 回分)。b = N^q は最下位 nibble に織り込み。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v4_7 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob24;
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
GNU_TARGET("vpclmulqdq") inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_shift= _mm256_srli_si256(_mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3)), 8);
 __m256i indices= _mm256_srli_si256(_mm256_srli_epi64(prod, 60), 8);
 return _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), _mm256_shuffle_epi8(RED_TABLE, indices));
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
 __m256i T78= mul2(T4, T34);
 tie(T[7], T[8])= unpack(T78);
 __m256i T910= mul2(T4, T56);
 tie(T[9], T[10])= unpack(T910);
 tie(T[11], T[12])= unpack(mul2(T4, T78));
 tie(T[13], T[14])= unpack(mul2(T4, T910));
 tie(T[15], B)= unpack(mul2(_mm256_set_epi64x(0, frob16(B), 0, T[7]), _mm256_set_epi64x(0, B, 0, T[8])));
 // メイン loop 2-lane 化: r = r_H·2^24 + r_L と分け a^r = frob24(a^{r_H})·a^{r_L} を
 // lockstep で計算 (直列チェーン init+11 反復 → init+5 反復 + frob24 + mul)。
 // ・chunk=0 の skip は廃止して無条件 mul2s (T[0]=1 の単位元掛け、分岐 mispredict も消える)
 // ・r ≥ 2^48 (確率 ~2^-16) のときだけ r_H が 7 nibble なので反復 +1 (r_L 側は先頭 0 詰め)
 // ・b = N^q は loop と独立に確定するので r_L の最下位 nibble operand に織り込み、
 //   末尾の直列を frob24 + mul 1 回に (bt は OoO が loop 実行中に計算)
 const u32 rl= r & 0xFFFFFF, rh= r >> 24;
 const u64 bt= mul(embed_idx(TABLES.PW_SIGMA_IDX[(u32(TABLES.LN_SIGMA[u16(B)]) * q) % 65535]), T[rl & 0xF]);
 const int it= int(r >> 48) + 5;
 u64 gl= T[(rl >> (4 * it)) & 0xF], gh= T[(rh >> (4 * it)) & 0xF];
 for(int i= it - 1; i >= 1; --i) tie(gh, gl)= unpack(mul2(_mm256_set_epi64x(0, frob4(gl), 0, frob4(gh)), _mm256_set_epi64x(0, T[(rl >> (4 * i)) & 0xF], 0, T[(rh >> (4 * i)) & 0xF])));
 tie(gh, gl)= unpack(mul2(_mm256_set_epi64x(0, frob4(gl), 0, frob4(gh)), _mm256_set_epi64x(0, bt, 0, T[rh & 0xF])));
 return mul(frob24(gh), gl);
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
