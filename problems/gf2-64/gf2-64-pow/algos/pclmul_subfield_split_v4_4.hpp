#pragma once
// pclmul_subfield_split_v4_2.hpp の preprocessing 一部 constexpr 化版:
//
// LN_SIGMA / PW_SIGMA_IDX (それぞれ 128 KiB) は σ^k chain (65535 muls) が必要で
// constexpr step limit を超えるので runtime init のまま。
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
// σ^0..σ^15
constexpr array<u64, 16> SIGMA_POW= {1ull, 11625825068197362226ull, 16726069499881173557ull, 1643712016162803251ull, 1139490258686435223ull, 2133510969792054338ull, 6902720057649047833ull, 6911235282440967732ull, 1609028687227055289ull, 17188600949917388119ull, 970962174813342580ull, 13487762679230420322ull, 17964985024964438593ull, 449778467375550729ull, 13547138401716955404ull, 5392541099451007413ull};
constexpr auto EMBED_BYTE= []() {
 array<array<u64, 256>, 2> embed_byte{};
 // contribution[k] = subfield 元 c で PEXT(c, MASK) = (1 << k) を満たすもの。
 // c = Σ c_i σ^i とおくと M·c_vec = e_k where M[r][i] = (σ^i >> picked[r]) & 1。
 // よって c_vec = M^{-1} · e_k = M^{-1} の k 列目。
 // 16x16 GF(2) 行列 [M | I] を Gauss-Jordan で [I | M^{-1}] に変形する。
 u32 M[16]= {}, Minv[16]= {};
 for(int r= 0; r < 16; ++r) {
  u32 row= 0;
  for(int i= 0; i < 16; ++i) {
   if((SIGMA_POW[i] >> r) & 1) row|= u32(1) << i;
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
 // contribution[k] = ⊕ over i: M^{-1}[i][k] == 1 of SIGMA_POW[i]
 // M^{-1}[i] (= Minv[i]) は 16-bit、その k bit が k 列目の i 行目成分
 u64 contribution[16]= {};
 for(int k= 0; k < 16; ++k) {
  u64 v= 0;
  for(int i= 0; i < 16; ++i) {
   if((Minv[i] >> k) & 1) v^= SIGMA_POW[i];
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
   embed_byte[p][b]= v;
  }
 }
 return embed_byte;
}();

// runtime-init される big tables (σ^k chain 経由)
inline u16 LN_SIGMA[65536];
inline u16 PW_SIGMA_IDX[65536];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 u64 cur= 1;
 for(u16 k= 0; k < 65535; ++k) {
  u16 idx= cur;
  LN_SIGMA[idx]= u16(k);
  PW_SIGMA_IDX[k]= u16(idx);
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[0]= 0;
 PW_SIGMA_IDX[65535]= 1;
}
inline u64 embed_idx(u16 idx) { return EMBED_BYTE[0][u8(idx)] ^ EMBED_BYTE[1][u8(idx >> 8)]; }
u64 pow_byte_window(u64 g, u64 e) {
 u64 T[16]= {1, g};
#pragma GCC unroll 14
 for(int i= 2; i < 16; ++i) T[i]= mul(T[i - 1], g);
 int top= 15 - (__builtin_clzll(e) >> 2);
 u64 acc= T[(e >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob4(acc);
  u32 chunk= u32((e >> (4 * i)) & 0xF);
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 if(!a) return 0;
 constexpr u64 M_VAL= (~u64(0)) / 65535u;
 const u16 q= e / M_VAL;
 if(!q) return pow_byte_window(a, e);
 const u64 r= e - M_VAL * q;
 const u16 N= mul(mul(a, frob16(a)), mul(frob32(a), frob48(a)));
 const u64 b= embed_idx(PW_SIGMA_IDX[(u32(LN_SIGMA[N]) * q) % 65535]);
 const u64 g= pow_byte_window(a, r);
 return mul(b, g);
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
