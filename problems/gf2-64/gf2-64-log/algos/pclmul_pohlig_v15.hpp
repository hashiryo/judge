#pragma once
// pclmul_pohlig_v7 + 共有 T precompute + VPCLMULQDQ:
//
// v7 の per-query は 4 つの pow_bw(x, EXP_*) を逐次呼び出し、 各々で T[16] precompute
// (14 mul) を行うため重複が大きい。 4 つとも base = x なので T は同一 → **1 回だけ** 計算。
//
// 加えて T precompute を VPCLMULQDQ の binary-tree で並列化:
//   L1: T[2]                                    (1 mul)
//   L2: T[3], T[4]                              (1 vpclmul)
//   L3: T[5..8]                                 (2 vpclmul)
//   L4: T[9..15]                                (3 vpclmul + 1 mul)
//   合計 ~8 ops (14 sequential mul から大幅短縮)
//
// 4 つの pow_bw_with_T は precomputed T を共有して main loop だけを回す。
//
// 効果見込み (per query):
//   元 v7: 4 × (T setup 14 mul + main loop ~10 mul) ≈ 4 × 24 = 96 mul
//   v8:    1 × T setup (8 ops via vpclmul) + 4 × main loop (~10 mul) ≈ 8 + 40 = 48 mul-equiv
//   ⇒ ~50% 短縮 per query。
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_log_pohlig_v11 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob8;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
// VPCLMUL 2 並列 mul + 並列 reduction (vmul_3_2 と同じ idiom)
VPCLMUL inline __m256i mul2(__m256i a_vec, __m256i b_vec, u64& r0, u64& r1) {
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
// T[i] = a^i for i = 0..15 を vpclmul binary-tree で構築
inline void build_T(u64 a, u64 T[16]) {
 T[0]= 1;
 T[1]= a;
 T[2]= mul(a, a);
 __m256i T2a= _mm256_set_epi64x(0, (i64)T[2], 0, (i64)a);
 // L2: T[3] = T[2]·a, T[4] = T[2]·T[2]
 __m256i T43= mul2(T2a, _mm256_set1_epi64x((i64)T[2]), T[3], T[4]);
 // L3: T[5..8]
 __m256i T4_v= _mm256_set1_epi64x((i64)T[4]);
 __m256i T56= mul2(T4_v, T2a, T[5], T[6]);
 mul2(T4_v, T43, T[7], T[8]);
 // L4: T[9..15]
 __m256i T8_v= _mm256_set1_epi64x((i64)T[8]);
 mul2(T8_v, T2a, T[9], T[10]);
 mul2(T8_v, T43, T[11], T[12]);
 mul2(T8_v, T56, T[13], T[14]);
 T[15]= mul(T[8], T[7]);
}
// precomputed T[16] を使った pow_bw (T 設定を skip するので main loop だけ)
inline u64 pow_bw_with_T(const u64 T[16], u64 e) {
 int top= 15;
 while(top > 0 && ((e >> (4 * top)) & 0xF) == 0) --top;
 u64 acc= T[(e >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob4(acc);
  u32 chunk= u32((e >> (4 * i)) & 0xF);
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
// init 用 (T 共有しない単発版)
inline u64 pow_bw(u64 a, u64 e) {
 u64 T[16];
 build_T(a, T);
 return pow_bw_with_T(T, e);
}
// =============================================================================
// F_{2^16}^* log table (compile-time 構築)
//
// 生成元 BETA = 0x1f1af3ec55a22e02 を natural 表現で chain (16-bit shift+XOR)。
// poly basis での低 16 bit が一意識別子になる選び方なので extract は単純な u16() cast。
// (vmul_3_2.hpp / subfield_split_constexpr_2.hpp と同じ idiom)
// 元 v8 の build_sigma_pext_table() (Gauss + 65535 muls runtime) を撤廃。
// =============================================================================
constexpr auto LN16= []() {
 // col[i] = u16(BETA^i) (poly basis 累乗の低 16 bit)
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
 array<u16, 65536> ln{};
 // BETA の最小多項式 y^16 + y^5 + y^3 + y^2 + 1, lower bits = 0x002D
 u16 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u16 lo= T_lo[u8(cur)] ^ T_hi[cur >> 8];
  ln[lo]= u32(u16(k)) * 2699 % 65535;  // H_LOG_INV を掛けておく (mul 1 回分の前処理)
  cur= u16(cur << 1) ^ (0x002DU & -u16(cur >> 15));
 }
 ln[0]= 0;
 return ln;
}();
// =============================================================================
// 65537 / 641 subgroup 用の direct hash log table
// =============================================================================
struct DirectLogTable {
 std::vector<u64> keys;
 std::vector<u32> values;
 u32 mask;
 void build(u64 base, u32 p) {
  u32 cap= 8;
  while(cap < 2u * p) cap*= 2;
  keys.assign(cap, ~u64(0));
  values.assign(cap, 0);
  mask= cap - 1;
  u64 cur= 1;
  for(u32 k= 0; k < p; ++k) {
   u32 h= cur & mask;
   while(keys[h] != ~u64(0)) h= (h + 1) & mask;
   keys[h]= cur;
   values[h]= k;
   cur= mul(cur, base);
  }
 }
 u32 lookup(u64 target) const {
  u32 h= target & mask;
  while(keys[h] != ~u64(0)) {
   if(keys[h] == target) return values[h];
   h= (h + 1) & mask;
  }
  return ~u32(0);
 }
};
inline DirectLogTable direct_641, direct_65537;
// =============================================================================
// BSGS for 6700417
// =============================================================================
struct BSGSTable {
 std::vector<u64> keys;
 std::vector<u32> values;
 u32 mask;
 u32 m;
 u64 q;
 u64 inv_base_m;
 void build(u64 base, u64 q_, u32 m_override= 0) {
  q= q_;
  if(m_override) m= m_override;
  else {
   m= 1;
   while((u64(m) * m) < q) ++m;
  }
  u32 cap= 8;
  while(cap < 4 * m) cap*= 2;
  keys.assign(cap, ~u64(0));
  values.assign(cap, 0);
  mask= cap - 1;
  u64 cur= 1;
  for(u64 j= 0; j < m; ++j) {
   u32 h= cur & mask;
   while(keys[h] != ~u64(0)) h= (h + 1) & mask;
   keys[h]= cur;
   values[h]= u32(j);
   cur= mul(cur, base);
  }
  inv_base_m= pow_bw(base, q - m);
 }
 u32 solve(u64 target) const {
  u64 t= target;
  u64 max_i= (q + m - 1) / m;
  for(u64 i= 0; i <= max_i; ++i) {
   u32 h= t & mask;
   while(keys[h] != ~u64(0)) {
    if(keys[h] == t) {
     u64 res= i * m + values[h];
     if(res < q) return u32(res);
    }
    h= (h + 1) & mask;
   }
   t= mul(t, inv_base_m);
  }
  return u32(q);
 }
};
inline BSGSTable bsgs_6700417;
constexpr u64 P_F16= 65535;
constexpr u64 P_641= 641;
constexpr u64 P_F17= 65537;
constexpr u64 P_BIG= 6700417;
constexpr u64 INV_65535_641= 243ULL;
constexpr u64 MOD2= (P_F16 * P_641) % P_F17;
constexpr u64 INV_MOD2_F17= 45242ULL;
constexpr u64 MOD3= (P_F16 * P_641 * P_F17) % P_BIG;
constexpr u64 INV_MOD3_BIG= 3883315ULL;
constexpr u64 MOD_F16= P_F16;
constexpr u64 MOD_F16_641= P_F16 * P_641;
constexpr u64 MOD_F16_641_F17= P_F16 * P_641 * P_F17;
constexpr u64 EXP_F16= 0x0001000100010001ull;
constexpr u64 EXP_641= 0x00663d80ff99c27full;
constexpr u64 EXP_F17= 0x0000ffff0000ffffull;
constexpr u64 EXP_BIG= 0x00000280fffffd7full;
// G_2^EXP_* は固定値なので Python (gf2_64.py) で事前計算して埋め込み:
//   g_641     = gf_pow(2, EXP_641)
//   g_65537   = gf_pow(2, EXP_F17)
//   g_6700417 = gf_pow(2, EXP_BIG)
constexpr u64 G_641= 0x6bf808f7824282a2ull;
constexpr u64 G_65537= 0x1c1e79669b95a7ceull;
constexpr u64 G_6700417= 0x00f542601703f991ull;
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;

 direct_641.build(G_641, u32(P_641));
 direct_65537.build(G_65537, u32(P_F17));
 bsgs_6700417.build(G_6700417, P_BIG, 131072);
}
u64 log_g(u64 x) {
 assert(x);
 u64 N= mul(x, frob32(x));
 const u16 x_f16= mul(N, frob16(N));
 u64 s= mul(x, sq(x));
 s= mul(s, frob2(s));
 s= mul(s, frob4(s));
 s= mul(s, frob8(s));  // 2^16-1
 const u64 x_65537= mul(s, frob32(s));
 u64 T[16];
 build_T(x, T);
 const u64 x_641= pow_bw_with_T(T, EXP_641);
 const u64 x_6700417= pow_bw_with_T(T, EXP_BIG);
 const u16 r1= LN16[x_f16];
 const u32 r0= direct_641.lookup(x_641);
 const u32 r2= direct_65537.lookup(x_65537);
 const u32 r3= bsgs_6700417.solve(x_6700417);
 const u16 cur_mod_641= r1 % P_641;
 const u16 diff0= (r0 + P_641 - cur_mod_641) % P_641;
 const u16 t0= (diff0 * INV_65535_641) % P_641;
 const u32 cur_after0_mod_F17= (r1 + MOD_F16 * t0) % P_F17;
 const u32 diff2= (r2 + P_F17 - cur_after0_mod_F17) % P_F17;
 const u32 t2= (diff2 * INV_MOD2_F17) % P_F17;
 const u32 cur_after2_mod_BIG= (r1 + (MOD_F16 * t0) % P_BIG + (MOD_F16_641 * t2) % P_BIG) % P_BIG;
 const u32 diff3= (r3 + P_BIG - cur_after2_mod_BIG) % P_BIG;
 const u32 t3= (u64(diff3) * INV_MOD3_BIG) % P_BIG;
 return u64(r1) + MOD_F16 * t0 + MOD_F16_641 * t2 + MOD_F16_641_F17 * t3;
}
}  // namespace gf2_64_log_pohlig_v11
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& xs) {
  using gf2_64_log_pohlig_v11::init_tables;
  using gf2_64_log_pohlig_v11::log_g;
  init_tables();
  vector<u64> ans(xs.size());
  for(size_t i= 0; i < xs.size(); ++i) ans[i]= log_g(xs[i]);
  return ans;
 }
};
