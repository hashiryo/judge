#pragma once
// frob_lookup と byte_window のハイブリッド:
//
//   α^e = ∏_{i=0..15, chunk_i ≠ 0} (T[chunk_i])^{2^{4i}}
//   T[c] = α^c  (c = 0..15, 14 mul で前計算)
//
// 各 (T[chunk_i])^{2^{4i}} は frob_{4i} byte table を 1 lookup するだけで得られる。
// 必要な byte table は k = 0, 4, 8, ..., 60 の 16 個 = 256 KiB (frob_lookup の 1 MiB より軽い)。
//
// byte_window との違い:
//   byte_window は frob4 を 15 回 serial に適用 → クリティカルパス長
//   本実装は 16 個の独立 frob_{4i} lookup → tree mul で depth log2 ≈ 4
// frob_lookup との違い:
//   frob_lookup は popcount(e) 個 (random で ~32) の lookup と mul
//   本実装は最大 16 個 (4-bit chunk 数) の lookup と mul → dense e で大幅に少ない
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_frob_lookup_byte_window {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
// FROB4K_BYTE[i][p][b] = (b << (8p))^{2^{4i}}, i = 0..15
inline u64 FROB4K_BYTE[16][8][256];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 for(int p= 0; p < 8; ++p)
  for(int b= 0; b < 256; ++b) FROB4K_BYTE[0][p][b]= u64(b) << (8 * p);
 for(int i= 1; i < 16; ++i)
  for(int p= 0; p < 8; ++p)
   for(int b= 0; b < 256; ++b) {
    u64 v= FROB4K_BYTE[i - 1][p][b];
    v= sq(v);
    v= sq(v);
    v= sq(v);
    v= sq(v);
    FROB4K_BYTE[i][p][b]= v;
   }
}
inline u64 apply_frob4k(int i, u64 a) {
 const auto& t= FROB4K_BYTE[i];
 return t[0][u8(a)] ^ t[1][u8(a >> 8)] ^ t[2][u8(a >> 16)] ^ t[3][u8(a >> 24)] ^ t[4][u8(a >> 32)] ^ t[5][u8(a >> 40)] ^ t[6][u8(a >> 48)] ^ t[7][u8(a >> 56)];
}
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 // T[c] = α^c, c = 0..15
 u64 T[16]= {1, a};
#pragma GCC unroll 14
 for(int j= 2; j < 16; ++j) T[j]= mul(T[j - 1], a);
 // 4-bit chunk ごとに frob_{4i}(T[chunk]) を収集
 u64 selected[16];
 int n= 0;
 for(int i= 0; i < 16; ++i) {
  u32 c= u32((e >> (4 * i)) & 0xF);
  if(c) selected[n++]= apply_frob4k(i, T[c]);
 }
 // Tree mul (depth ≤ 4)
 while(n > 1) {
  int new_n= 0, j= 0;
  for(; j + 1 < n; j+= 2) selected[new_n++]= mul(selected[j], selected[j + 1]);
  if(j < n) selected[new_n++]= selected[j];
  n= new_n;
 }
 return selected[0];
}
}  // namespace gf2_64_pow_frob_lookup_byte_window
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_frob_lookup_byte_window::init_tables;
  using gf2_64_pow_frob_lookup_byte_window::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
