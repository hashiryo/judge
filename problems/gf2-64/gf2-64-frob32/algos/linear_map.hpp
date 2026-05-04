#pragma once
// a^32 を F_2-線型写像と見て、 8 byte tables (16 KB) の lookup で計算。
//   a^32 = XOR_{p=0..7} FROB_BYTE[p][(a >> 8p) & 0xFF]
//
// init: 8 × 256 = 2048 個の sq^5(b << 8p) を runtime 計算。
// per query: 8 lookup + 8 XOR ≈ 8-12 cyc。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
namespace gf2_64_frob32_linear_map {
inline u64 sq_naive(u64 a) {
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
inline u64 frob32_naive(u64 a) {
 for(int i= 0; i < 5; ++i) a= sq_naive(a);
 return a;
}
inline u64 FROB_BYTE[8][256];
inline bool inited= false;
inline void init_tables() {
 if(inited) return;
 inited= true;
 for(int p= 0; p < 8; ++p)
  for(int b= 0; b < 256; ++b) FROB_BYTE[p][b]= frob32_naive(u64(b) << (8 * p));
}
inline u64 frob32(u64 a) {
 return FROB_BYTE[0][u8(a)] ^ FROB_BYTE[1][u8(a >> 8)] ^ FROB_BYTE[2][u8(a >> 16)] ^ FROB_BYTE[3][u8(a >> 24)] ^ FROB_BYTE[4][u8(a >> 32)] ^ FROB_BYTE[5][u8(a >> 40)] ^ FROB_BYTE[6][u8(a >> 48)] ^ FROB_BYTE[7][u8(a >> 56)];
}
}

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using namespace gf2_64_frob32_linear_map;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= frob32(as[i]);
  return ans;
 }
};
