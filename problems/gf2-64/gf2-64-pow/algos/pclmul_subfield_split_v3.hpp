#pragma once
// pclmul_subfield_split_v2.hpp の発展形:
// F_{2^16}^* の log/exp を Nimber 基底 (basis_change の poly_to_nim) から
// σ-coeff 基底 + PEXT に置換。
//
//   Extract: PEXT(N, MASK)        — 1 命令
//   Embed:   PW_SIGMA[log]        — 直接 64-bit poly を返すテーブル (512 KiB)
//   LN:      LN_SIGMA[pext_idx]   — 16-bit log (128 KiB)
//
// メモリ: 640 KiB (L2 borderline、spill 時は L3 hit) - v2 より重い
// Per-pow: F_{2^16} 操作が 8 byte lookup → PEXT 1 命令で短縮
//
// Nimber 依存ゼロ。div の pclmul_norm_pext.hpp と同じ基底思想。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#define HAVE_PEXT 1
#else
#define PCLMUL_RUN
#define HAVE_PEXT 0
#endif
namespace gf2_64_pow_subfield_split_v3 {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob48;
// =============================================================================
// σ-coeff 基底 + PEXT
// =============================================================================
constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;

inline u64 PEXT_MASK= 0;
inline u16 LN_SIGMA[65536];  // PEXT idx → log (16-bit, 128 KiB)
inline u64 PW_SIGMA[65535];  // log → 64-bit poly (= σ^log) directly (512 KiB)

#if !HAVE_PEXT
inline int PEXT_POS[16];
#endif
[[gnu::target("pclmul")]] void build_sigma_tables() {
 // σ^0..σ^15 を計算 → 線型独立な 16 bit 位置を Gauss 消去で見つける
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
#if !HAVE_PEXT
 for(int i= 0; i < 16; ++i) PEXT_POS[i]= picked[i];
#endif
 // LN_SIGMA, PW_SIGMA を σ^k chain で構築
 u64 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u32 idx;
#if HAVE_PEXT
  idx= u32(_pext_u64(cur, PEXT_MASK));
#else
  idx= 0;
  for(int i= 0; i < 16; ++i) idx|= u32((cur >> picked[i]) & 1) << i;
#endif
  LN_SIGMA[idx]= u16(k);
  PW_SIGMA[k]= cur;
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[0]= 0;  // for x=0 case (shouldn't happen)
}
[[gnu::always_inline]] inline u32 extract_idx(u64 N) {
#if HAVE_PEXT
 return u32(_pext_u64(N, PEXT_MASK));
#else
 u32 r= 0;
 for(int i= 0; i < 16; ++i) r|= u32((N >> PEXT_POS[i]) & 1) << i;
 return r;
#endif
}
[[gnu::always_inline]] inline u32 e_mod_65535(u64 e) {
 const u32 s= u32(e & 0xFFFF) + u32((e >> 16) & 0xFFFF) + u32((e >> 32) & 0xFFFF) + u32((e >> 48) & 0xFFFF);
 u32 r= s;
 if(r >= 65535) r-= 65535;
 if(r >= 65535) r-= 65535;
 if(r >= 65535) r-= 65535;
 return r;
}
[[gnu::target("pclmul")]] u64 pow_byte_window(u64 g, u64 e) {
 if(e == 0) return 1;
 u64 T[16];
 T[0]= 1;
 T[1]= g;
#pragma GCC unroll 14
 for(int i= 2; i < 16; ++i) T[i]= mul(T[i - 1], g);
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
constexpr u32 M_INV_MOD_65535= 16384;
inline bool inited= false;
[[gnu::target("pclmul")]] void init_tables() {
 if(inited) return;
 inited= true;
 // frob byte tables は _shared/frob.hpp が compile-time に提供
 build_sigma_tables();
}
[[gnu::target("pclmul")]] u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 // ① N(α) 計算 (3 frob 並列)
 const u64 a16= frob16(a);
 const u64 a32= frob32(a);
 const u64 a48= frob48(a);
 const u64 N= mul(mul(a, a16), mul(a32, a48));
 // ② β 抽出 via PEXT + log/exp
 const u32 N_idx= extract_idx(N);
 if(N_idx == 0) return 0;  // α = 0 (or N=0)
 const u32 log_N= LN_SIGMA[N_idx];
 const u32 log_beta= u32((u64(log_N) * M_INV_MOD_65535) % 65535);
 // ③ γ = α · β^{-1} via PW_SIGMA 直接 lookup
 const u32 log_beta_inv= (65535u - log_beta) % 65535u;
 const u64 beta_inv_poly= PW_SIGMA[log_beta_inv];
 const u64 gamma= mul(a, beta_inv_poly);
 // ④ e の分解
 constexpr u64 M_VAL= (~u64(0)) / 65535u;
 const u32 e_low_red= e_mod_65535(e);
 const u64 e_high= e % M_VAL;
 // ⑤ β^{e_low}
 const u64 beta_pow= (e_low_red == 0) ? 1ull : PW_SIGMA[u32((u64(log_beta) * e_low_red) % 65535)];
 // ⑥ γ^{e_high}
 const u64 gamma_pow= pow_byte_window(gamma, e_high);
 // ⑦ 結合
 return mul(beta_pow, gamma_pow);
}
}  // namespace gf2_64_pow_subfield_split_v3
struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v3::init_tables;
  using gf2_64_pow_subfield_split_v3::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
