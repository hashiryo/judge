#pragma once

#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pclmul_itoh_tsujii {
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
u64 inv(u64 a) {
 constexpr u64 e= 0xFFFFFFFFFFFFFFFEull;  // 2^64 - 2
 // Precompute T[i] = a^i for i = 0..15 (14 muls)
 u64 T[16]= {1, a};
 for(int i= 2; i < 16; ++i) T[i]= mul(T[i - 1], a);

 u64 acc= T[15];
 for(int i= 14; i >= 0; --i) {
  acc= frob4(acc);
  u16 chunk= (e >> (4 * i)) & 0xF;
  acc= mul(acc, T[chunk]);
 }
 return acc;
}
}  // namespace gf2_64_pclmul_itoh_tsujii
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_itoh_tsujii::inv;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
