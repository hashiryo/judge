#pragma once

#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pclmul_window {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// constexpr で [2][256] テーブルを構築
constexpr auto EMBED_BYTE= []() {
 constexpr std::array<uint64_t, 16> SUBFIELD_BASIS= {0x0000000000000001ULL, 0x5fbfaec6aeac0002ULL, 0xb06c601895640004ULL, 0xb013b5277b7c0008ULL, 0xb5ebb915248a0010ULL, 0x109bb25b2c600020ULL, 0xbf3bd95bd4190040ULL, 0x0fc66342279b0080ULL, 0xb6418f5e57c50100ULL, 0xaa194bd4b83f0200ULL, 0x1b5217b4dcc70400ULL, 0xbb06fa73867a0800ULL, 0x006fd55b23331000ULL, 0x4ae8fb39198c2000ULL, 0xfbd141b29b4f4000ULL, 0x1d9ce1776be78000ULL};
 std::array<std::array<uint64_t, 256>, 2> T{};
 for(int half= 0; half < 2; ++half) {
  for(int i= 0; i < 256; ++i) {
   uint64_t v= 0;
   for(int b= 0; b < 8; ++b) {
    if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
   }
   T[half][i]= v;
  }
 }
 return T;
}();
inline u64 embed_idx(u16 idx) { return EMBED_BYTE[0][u8(idx)] ^ EMBED_BYTE[1][u8(idx >> 8)]; }
// runtime-init される big tables (σ^k chain 経由)
inline u16 INV[65535];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 u16 LN_SIGMA[65536];
 u16 PW_SIGMA_IDX[65535];
 constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;
 u64 cur= 1;
 for(u16 k= 0; k < 65535; ++k) {
  u16 idx= cur;
  LN_SIGMA[idx]= 65535 - k;
  PW_SIGMA_IDX[k]= idx;
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[1]= 0;
 for(u16 k= 1; k < 65535; ++k) INV[k]= PW_SIGMA_IDX[LN_SIGMA[k]];
}
u64 inv(u64 a) {
 assert(a != 0);
 u64 a16= frob16(a);
 u64 a32= frob16(a16);
 u64 a48= frob16(a32);
 u64 g= mul(a16, mul(a32, a48));
 u64 b= embed_idx(INV[u16(mul(g, a))]);
 return mul(b, g);
}
}  // namespace gf2_64_pclmul_window
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_window::init_tables;
  using gf2_64_pclmul_window::inv;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
