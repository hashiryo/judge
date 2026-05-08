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
namespace gf2_64_log_pohlig_v8 {
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
// VPCLMUL 2 並列 mul + 並列 reduction (vmul_3_2 と同じ idiom)
inline __m256i mul2(__m256i a_vec, __m256i b_vec, u64& r0, u64& r1) {
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
 if(e == 0) return 1;
 int top= 15;
 while(top > 0 && ((e >> (4 * top)) & 0xF) == 0) --top;
 u64 acc= T[(e >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob4(acc);
  unsigned chunk= unsigned((e >> (4 * i)) & 0xF);
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
// init 用 (T 共有しない単発版)
inline u64 pow_bw(u64 a, u64 e) {
 if(e == 0) return 1;
 u64 T[16];
 build_T(a, T);
 return pow_bw_with_T(T, e);
}
// =============================================================================
// σ + PEXT for F_{2^16}^* (v7 と同じ)
// =============================================================================
constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;
inline u64 PEXT_MASK= 0;
inline u16 LN_SIGMA_PEXT[65536];
void build_sigma_pext_table() {
 u64 sigma_pow[16];
 sigma_pow[0]= 1;
 for(int i= 1; i < 16; ++i) sigma_pow[i]= mul(sigma_pow[i - 1], SIGMA);
 u16 row_vec[64];
 for(int r= 0; r < 64; ++r) {
  u16 v= 0;
  for(int c= 0; c < 16; ++c) {
   if((sigma_pow[c] >> r) & 1) v|= u16(1) << c;
  }
  row_vec[r]= v;
 }
 int picked[16];
 int n_picked= 0;
 u16 basis[16]= {};
 for(int r= 0; r < 64 && n_picked < 16; ++r) {
  u16 v= row_vec[r];
  for(int k= 15; k >= 0 && v; --k) {
   if(!((v >> k) & 1)) continue;
   if(basis[k] == 0) {
    basis[k]= v;
    picked[n_picked++]= r;
    break;
   }
   v^= basis[k];
  }
 }
 PEXT_MASK= 0;
 for(int i= 0; i < 16; ++i) PEXT_MASK|= (u64(1) << picked[i]);
 u64 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u32 idx= u32(_pext_u64(cur, PEXT_MASK));
  LN_SIGMA_PEXT[idx]= u16(k);
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA_PEXT[0]= 0;
}
inline u32 extract_idx(u64 N) { return u32(_pext_u64(N, PEXT_MASK)); }
// =============================================================================
// 65537 / 641 subgroup 用の direct hash log table
// =============================================================================
struct DirectLogTable {
 std::vector<u64> keys;
 std::vector<u32> values;
 u64 mask;
 void build(u64 base, u32 p) {
  u64 cap= 8;
  while(cap < 2u * p) cap*= 2;
  keys.assign(cap, ~u64(0));
  values.assign(cap, 0);
  mask= cap - 1;
  u64 cur= 1;
  for(u32 k= 0; k < p; ++k) {
   u64 h= (cur * 0x9E3779B97F4A7C15ull) & mask;
   while(keys[h] != ~u64(0)) h= (h + 1) & mask;
   keys[h]= cur;
   values[h]= k;
   cur= mul(cur, base);
  }
 }
 u32 lookup(u64 target) const {
  u64 h= (target * 0x9E3779B97F4A7C15ull) & mask;
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
 u64 mask;
 u64 m;
 u64 q;
 u64 inv_base_m;
 void build(u64 base, u64 q_, u64 m_override= 0) {
  q= q_;
  if(m_override) m= m_override;
  else {
   m= 1;
   while(m * m < q) ++m;
  }
  u64 cap= 8;
  while(cap < 4 * m) cap*= 2;
  keys.assign(cap, ~u64(0));
  values.assign(cap, 0);
  mask= cap - 1;
  u64 cur= 1;
  for(u64 j= 0; j < m; ++j) {
   u64 h= (cur * 0x9E3779B97F4A7C15ull) & mask;
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
   u64 h= (t * 0x9E3779B97F4A7C15ull) & mask;
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
inline u32 H_LOG_INV;
u32 solve_f16(u64 x_proj_poly) {
 const u32 idx= extract_idx(x_proj_poly);
 if(idx == 0) return 0;
 const u32 log_x= LN_SIGMA_PEXT[idx];
 return u32((u64(log_x) * H_LOG_INV) % 65535);
}
constexpr u64 modinv_runtime(u64 a, u64 m) {
 long long r0= (long long)m, r1= (long long)a;
 long long s0= 0, s1= 1;
 while(r1 != 0) {
  long long q= r0 / r1;
  long long r2= r0 - q * r1;
  r0= r1;
  r1= r2;
  long long s2= s0 - q * s1;
  s0= s1;
  s1= s2;
 }
 if(s0 < 0) s0+= (long long)m;
 return (u64)s0;
}
constexpr u64 P_F16= 65535;
constexpr u64 P_641= 641;
constexpr u64 P_F17= 65537;
constexpr u64 P_BIG= 6700417;
constexpr u64 INV_65535_641= modinv_runtime(65535 % 641, 641);
constexpr u64 MOD2= (P_F16 * P_641) % P_F17;
constexpr u64 INV_MOD2_F17= modinv_runtime(MOD2, P_F17);
constexpr u64 MOD3= (P_F16 * P_641 * P_F17) % P_BIG;
constexpr u64 INV_MOD3_BIG= modinv_runtime(MOD3, P_BIG);
constexpr u64 MOD_F16= P_F16;
constexpr u64 MOD_F16_641= P_F16 * P_641;
constexpr u64 MOD_F16_641_F17= P_F16 * P_641 * P_F17;
constexpr u64 G_2= 2;
constexpr u64 EXP_F16= 0x0001000100010001ull;
constexpr u64 EXP_641= 0x00663d80ff99c27full;
constexpr u64 EXP_F17= 0x0000ffff0000ffffull;
constexpr u64 EXP_BIG= 0x00000280fffffd7full;
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 build_sigma_pext_table();
 const u64 h_poly= pow_bw(G_2, EXP_F16);
 const u32 h_idx= extract_idx(h_poly);
 const u32 h_log= LN_SIGMA_PEXT[h_idx];
 H_LOG_INV= u32(modinv_runtime(h_log, 65535));
 const u64 g_641= pow_bw(G_2, EXP_641);
 const u64 g_65537= pow_bw(G_2, EXP_F17);
 const u64 g_6700417= pow_bw(G_2, EXP_BIG);
 direct_641.build(g_641, u32(P_641));
 direct_65537.build(g_65537, u32(P_F17));
 bsgs_6700417.build(g_6700417, P_BIG, 131072);
}
u64 log_g(u64 x) {
 // T[16] を 1 回だけ vpclmul tree で構築 → 4 つの pow_bw で共有
 u64 T[16];
 build_T(x, T);
 const u64 x_f16= pow_bw_with_T(T, EXP_F16);
 const u64 x_641= pow_bw_with_T(T, EXP_641);
 const u64 x_65537= pow_bw_with_T(T, EXP_F17);
 const u64 x_6700417= pow_bw_with_T(T, EXP_BIG);
 const u32 r1= solve_f16(x_f16);
 const u32 r0= direct_641.lookup(x_641);
 const u32 r2= direct_65537.lookup(x_65537);
 const u32 r3= bsgs_6700417.solve(x_6700417);
 const u64 cur_mod_641= u64(r1) % P_641;
 const u64 diff0= (r0 + P_641 - cur_mod_641) % P_641;
 const u64 t0= (diff0 * INV_65535_641) % P_641;
 const u64 cur_after0_mod_F17= (u64(r1) + MOD_F16 * t0) % P_F17;
 const u64 diff2= (r2 + P_F17 - cur_after0_mod_F17) % P_F17;
 const u64 t2= (diff2 * INV_MOD2_F17) % P_F17;
 const u64 cur_after2_mod_BIG= (u64(r1) + (MOD_F16 * t0) % P_BIG + (MOD_F16_641 * t2) % P_BIG) % P_BIG;
 const u64 diff3= (r3 + P_BIG - cur_after2_mod_BIG) % P_BIG;
 const u64 t3= (__uint128_t(diff3) * INV_MOD3_BIG) % P_BIG;
 return u64(r1) + MOD_F16 * t0 + MOD_F16_641 * t2 + MOD_F16_641_F17 * t3;
}
}  // namespace gf2_64_log_pohlig_v8
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& xs) {
  using gf2_64_log_pohlig_v8::init_tables;
  using gf2_64_log_pohlig_v8::log_g;
  init_tables();
  vector<u64> ans(xs.size());
  for(size_t i= 0; i < xs.size(); ++i) ans[i]= log_g(xs[i]);
  return ans;
 }
};
