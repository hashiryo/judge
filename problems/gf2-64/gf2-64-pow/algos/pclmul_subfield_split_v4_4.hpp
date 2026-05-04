#pragma once
// pclmul_subfield_split_v4_2.hpp の preprocessing 一部 constexpr 化版:
//
// LN_SIGMA / PW_SIGMA_IDX (それぞれ 128 KiB) は σ^k chain (65535 muls) が必要で
// constexpr step limit を超えるので runtime init のまま。
// PEXT_MASK と EMBED_BYTE (4 KiB) は σ^0..σ^15 (16 muls) と Gauss/逆行列だけで
// 決まるので compile-time に押し込む。
//
// 利点: 細かい table 構築コードと frob byte chain の組み合わせ build を分離、init
//   コストは σ^k chain の純コストだけになる (元 v4_2: chain + EMBED 構築 + Gauss)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v4_3 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob48;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;

namespace _ce {
// constexpr GF(2^64) mul: P(x) = x^64 + x^4 + x^3 + x + 1, RED = 0x1B
constexpr u64 mul_gf264(u64 a, u64 b) {
 u64 lo= 0, hi= 0;
 while(b) {
  int i= __builtin_ctzll(b);
  b&= b - 1;
  lo^= a << i;
  if(i) hi^= a >> (64 - i);
 }
 while(hi) {
  u64 lo2= 0, hi2= 0;
  u64 h= hi;
  while(h) {
   int i= __builtin_ctzll(h);
   h&= h - 1;
   lo2^= 0x1Bull << i;
   if(i >= 60) hi2^= 0x1Bull >> (64 - i);
  }
  lo^= lo2;
  hi= hi2;
 }
 return lo;
}
struct StaticTables {
 u64 PEXT_MASK;
 int picked[16];
 u64 EMBED_BYTE[2][256];  // 16-bit subfield idx → 64-bit poly
};
constexpr StaticTables make_static_tables() {
 StaticTables t{};
 // σ^0..σ^15
 u64 sigma_pow[16]= {};
 sigma_pow[0]= 1;
 for(int i= 1; i < 16; ++i) sigma_pow[i]= mul_gf264(sigma_pow[i - 1], SIGMA);
 // Gauss 消去で線型独立 16 bit 位置を特定
 u16 row_vec[64]= {};
 for(int r= 0; r < 64; ++r) {
  u16 v= 0;
  for(int c= 0; c < 16; ++c) {
   if((sigma_pow[c] >> r) & 1) v|= u16(1) << c;
  }
  row_vec[r]= v;
 }
 int n_picked= 0;
 u16 basis[16]= {};
 for(int r= 0; r < 64 && n_picked < 16; ++r) {
  u16 v= row_vec[r];
  for(int k= 15; k >= 0 && v; --k) {
   if(!((v >> k) & 1)) continue;
   if(basis[k] == 0) {
    basis[k]= v;
    t.picked[n_picked++]= r;
    break;
   }
   v^= basis[k];
  }
 }
 t.PEXT_MASK= 0;
 for(int i= 0; i < 16; ++i) t.PEXT_MASK|= (u64(1) << t.picked[i]);
 // contribution[k] = subfield 元 c で PEXT(c, MASK) = (1 << k) を満たすもの。
 // c = Σ c_i σ^i とおくと M·c_vec = e_k where M[r][i] = (σ^i >> picked[r]) & 1。
 // よって c_vec = M^{-1} · e_k = M^{-1} の k 列目。
 // 16x16 GF(2) 行列 [M | I] を Gauss-Jordan で [I | M^{-1}] に変形する。
 u32 M[16]= {}, Minv[16]= {};
 for(int r= 0; r < 16; ++r) {
  u32 row= 0;
  for(int i= 0; i < 16; ++i) {
   if((sigma_pow[i] >> t.picked[r]) & 1) row|= u32(1) << i;
  }
  M[r]= row;
  Minv[r]= u32(1) << r;
 }
 for(int col= 0; col < 16; ++col) {
  int pv= -1;
  for(int r= col; r < 16; ++r)
   if((M[r] >> col) & 1) {
    pv= r;
    break;
   }
  if(pv == -1) continue;
  if(pv != col) {
   u32 tm= M[col];
   M[col]= M[pv];
   M[pv]= tm;
   tm= Minv[col];
   Minv[col]= Minv[pv];
   Minv[pv]= tm;
  }
  for(int r= 0; r < 16; ++r) {
   if(r != col && ((M[r] >> col) & 1)) {
    M[r]^= M[col];
    Minv[r]^= Minv[col];
   }
  }
 }
 // contribution[k] = ⊕ over i: M^{-1}[i][k] == 1 of sigma_pow[i]
 // M^{-1}[i] (= Minv[i]) は 16-bit、その k bit が k 列目の i 行目成分
 u64 contribution[16]= {};
 for(int k= 0; k < 16; ++k) {
  u64 v= 0;
  for(int i= 0; i < 16; ++i) {
   if((Minv[i] >> k) & 1) v^= sigma_pow[i];
  }
  contribution[k]= v;
 }
 // EMBED_BYTE: 16-bit idx の low/high byte に対する寄与
 for(int p= 0; p < 2; ++p) {
  for(int b= 0; b < 256; ++b) {
   u64 v= 0;
   for(int bit= 0; bit < 8; ++bit) {
    if((b >> bit) & 1) v^= contribution[p * 8 + bit];
   }
   t.EMBED_BYTE[p][b]= v;
  }
 }
 return t;
}
}  // namespace _ce
inline constexpr auto STATIC_TABLES= _ce::make_static_tables();

// runtime-init される big tables (σ^k chain 経由)
inline u16 LN_SIGMA[65536];
inline u16 PW_SIGMA_IDX[65535];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 u64 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u32 idx= u32(_pext_u64(cur, STATIC_TABLES.PEXT_MASK));
  LN_SIGMA[idx]= u16(k);
  PW_SIGMA_IDX[k]= u16(idx);
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[0]= 0;
}

inline u32 extract_idx(u64 N) { return u32(_pext_u64(N, STATIC_TABLES.PEXT_MASK)); }
inline u64 embed_idx(u16 idx) { return STATIC_TABLES.EMBED_BYTE[0][u8(idx)] ^ STATIC_TABLES.EMBED_BYTE[1][u8(idx >> 8)]; }
inline u32 e_mod_65535(u64 e) {
 u32 s= u16(e) + u16(e >> 16) + u16(e >> 32) + u16(e >> 48);
 while(s >= 65535) s-= 65535;
 return s;
}
u64 pow_byte_window(u64 g, u64 e) {
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
u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 const u64 N= mul(mul(a, frob16(a)), mul(frob32(a), frob48(a)));
 const u32 N_idx= extract_idx(N);
 if(N_idx == 0) return 0;
 constexpr u64 M_VAL= (~u64(0)) / 65535u;
 const u32 log_N= LN_SIGMA[N_idx];
 const u32 log_beta= u32((u64(log_N) * M_INV_MOD_65535) % 65535);
 const u32 log_beta_inv= (65535u - log_beta) % 65535u;
 const u64 beta_inv_poly= embed_idx(PW_SIGMA_IDX[log_beta_inv]);
 const u64 gamma= mul(a, beta_inv_poly);
 const u32 e_low_red= e_mod_65535(e);
 const u64 e_high= e % M_VAL;
 const u64 beta_pow= (e_low_red == 0) ? 1ull : embed_idx(PW_SIGMA_IDX[u32((u64(log_beta) * e_low_red) % 65535)]);
 const u64 gamma_pow= pow_byte_window(gamma, e_high);
 return mul(beta_pow, gamma_pow);
}
}  // namespace gf2_64_pow_subfield_split_v4_3
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v4_3::init_tables;
  using gf2_64_pow_subfield_split_v4_3::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
