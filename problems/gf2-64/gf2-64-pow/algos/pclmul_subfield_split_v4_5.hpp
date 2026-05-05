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
inline u16 LN_SIGMA[65536];
inline u16 PW_SIGMA_IDX[65535];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;
 u64 cur= 1;
 for(u16 k= 0; k < 65535; ++k) {
  u16 idx= cur;
  LN_SIGMA[idx]= k;
  PW_SIGMA_IDX[k]= u16(idx);
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[0]= 0;
}
u64 pow_byte_window(u64 g, u64 e) {
 u64 T[16]= {1, g};
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
 const u16 N= mul(mul(a, frob16(a)), mul(frob32(a), frob48(a)));
 const u64 b= embed_idx(PW_SIGMA_IDX[(u32(LN_SIGMA[N]) * q) % 65535]);
 const u64 r= e - M_VAL * q;
 if(!r) return b;
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
