#pragma once
// a^16 を F_2-線型写像と見て、 8 byte tables (16 KB) の lookup で計算。
//   a^16 = XOR_{p=0..7} FROB_BYTE[p][(a >> 8p) & 0xFF]
//
// **constexpr 化**: PDEP intrinsic は constexpr 不可だが、 等価な bit-interleave
// shift + AND チェインなら constexpr 評価できる。 sq_naive (bit-by-bit、 ~512 ops/call)
// から fast sq (~25 ops/call) に置き換えて step 上限内に収める。
//
// per query: 8 lookup + 8 XOR ≈ 8-12 cyc。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
namespace gf2_64_frob16_linear_map {
constexpr array<array<u64, 256>, 8> FROB_BYTE= [] {
 // PDEP fallback と同じ bit interleave (constexpr-ok)
 auto spread= [](u32 a) constexpr -> u64 {
  u64 x= a;
  x= (x | (x << 16)) & 0x0000FFFF0000FFFFull;
  x= (x | (x << 8)) & 0x00FF00FF00FF00FFull;
  x= (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
  x= (x | (x << 2)) & 0x3333333333333333ull;
  x= (x | (x << 1)) & 0x5555555555555555ull;
  return x;
 };
 constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
 auto sq= [&](u64 a) constexpr -> u64 {
  u64 h= spread(u32(a >> 32));
  u64 d= h ^ (h << 1);
  return spread(u32(a)) ^ RED[h >> 60] ^ d ^ (d << 3);
 };
 array<array<u64, 256>, 8> t{};
 for(int p= 0; p < 8; ++p) {
  for(int b= 0; b < 256; ++b) {
   u64 v= u64(b) << (8 * p);
   for(int i= 0; i < 4; ++i) v= sq(v);
   t[p][b]= v;
  }
 }
 return t;
}();
inline u64 frob16(u64 a) {
 return FROB_BYTE[0][u8(a)] ^ FROB_BYTE[1][u8(a >> 8)] ^ FROB_BYTE[2][u8(a >> 16)] ^ FROB_BYTE[3][u8(a >> 24)] ^ FROB_BYTE[4][u8(a >> 32)] ^ FROB_BYTE[5][u8(a >> 40)] ^ FROB_BYTE[6][u8(a >> 48)] ^ FROB_BYTE[7][u8(a >> 56)];
}
}

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_frob16_linear_map::frob16;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= frob16(as[i]);
  return ans;
 }
};
