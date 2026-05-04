#pragma once
// Frobenius word-table linear map:
//   入力を 4 個の 16-bit chunk に分けて、 各位置 p ∈ {0,1,2,3} ごとに 65536-entry の
//   lookup table を precompute。
//   sq(a) = XOR_{p=0..3} FROB1_WORD[p][(a >> 16p) & 0xFFFF]
//
//   per query: 4 lookups + 4 XOR ≈ 4-6 cyc (byte 版の半分)
//
// trade-off vs byte 版 ([8][256] = 16KB):
//   - lookup 数: 8 → 4 (半減)
//   - メモリ: 16 KB (L1 fit) → 2 MB (L1/L2 spill, L3 hit)
//   - random access pattern なので cache miss 増加リスク。 L3 hit ~30 cyc / lookup
//   - per-query: byte 版の 8 lookups × 5 cyc = 40 cyc vs word 版の 4 lookups × 30 cyc = 120 cyc
//     ↑ word が遅くなる可能性も。 実測で勝負
//
// constexpr で compile-time 構築 (init unnecessary)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
namespace gf2_64_sq_frobenius_word {
constexpr u64 sq_naive(u64 a) {
 u64 lo= 0, hi= 0;
 for(int i= 0; i < 64; ++i) {
  if((a >> i) & 1) {
   int j= 2 * i;
   if(j < 64) lo^= u64(1) << j;
   else hi^= u64(1) << (j - 64);
  }
 }
 for(int i= 63; i >= 0; --i) {
  if((hi >> i) & 1) {
   hi^= u64(1) << i;
   lo^= IRRED_LOW << i;
   if(i > 0) hi^= IRRED_LOW >> (64 - i);
  }
 }
 return lo;
}
inline const auto FROB1_WORD= [] {
 // 2 MB. constexpr で持つと compile 時間爆発するので runtime init。
 // (compile time eval なら byte 版同様の constexpr 化でも良い)
 static array<array<u64, 65536>, 4> t{};
 for(int p= 0; p < 4; ++p)
  for(int w= 0; w < 65536; ++w) t[p][w]= sq_naive(u64(w) << (16 * p));
 return &t;
}();
inline u64 sq(u64 a) {
 return (*FROB1_WORD)[0][u16(a)] ^ (*FROB1_WORD)[1][u16(a >> 16)] ^ (*FROB1_WORD)[2][u16(a >> 32)] ^ (*FROB1_WORD)[3][u16(a >> 48)];
}
}
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_frobenius_word::sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sq(as[i]);
  return ans;
 }
};
