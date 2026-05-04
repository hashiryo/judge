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
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_frob_lookup_byte_window {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
namespace _detail {
struct FrobBytes {
 u64 t[16][8][256];
};
// 既存 constexpr FROB4_BYTE (std::array) を 1 度だけ raw 配列に転写してから使うと
// constexpr step が抑えられる (libc++ の array access overhead 回避)
constexpr FrobBytes make_frob4k_byte() {
 FrobBytes r{};
 // ローカル raw 配列に FROB4_BYTE を転写 (以降 raw access のみで step 節約)
 u64 F4[8][256]{};
 for(int p= 0; p < 8; ++p)
  for(int b= 0; b < 256; ++b) F4[p][b]= gf2_64_pclmul::FROB4_BYTE[p][b];
 // i = 0: identity
 for(int p= 0; p < 8; ++p)
  for(int b= 0; b < 256; ++b) r.t[0][p][b]= u64(b) << (8 * p);
 // i > 0: r.t[i] = frob4 ∘ r.t[i-1]  (raw 配列で apply_frob4)
 for(int i= 1; i < 16; ++i)
  for(int p= 0; p < 8; ++p)
   for(int b= 0; b < 256; ++b) {
    u64 a= r.t[i - 1][p][b];
    r.t[i][p][b]= F4[0][u8(a)] ^ F4[1][u8(a >> 8)] ^ F4[2][u8(a >> 16)] ^ F4[3][u8(a >> 24)] ^ F4[4][u8(a >> 32)] ^ F4[5][u8(a >> 40)] ^ F4[6][u8(a >> 48)] ^ F4[7][u8(a >> 56)];
   }
 return r;
}
}  // namespace _detail
inline constexpr auto FROB4K_BYTE_S= _detail::make_frob4k_byte();
inline constexpr auto& FROB4K_BYTE= FROB4K_BYTE_S.t;

inline u64 apply_frob4k(int i, u64 a) {
 const auto& t= FROB4K_BYTE[i];
 return t[0][u8(a)] ^ t[1][u8(a >> 8)] ^ t[2][u8(a >> 16)] ^ t[3][u8(a >> 24)] ^ t[4][u8(a >> 32)] ^ t[5][u8(a >> 40)] ^ t[6][u8(a >> 48)] ^ t[7][u8(a >> 56)];
}
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
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
  using gf2_64_pow_frob_lookup_byte_window::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
