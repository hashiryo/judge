#pragma once

#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pclmul_window {
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
u64 inv(u64 a) {
 // Precompute T[i] = a^i for i = 0..15 (14 muls)
 u64 T14= sq(mul(mul(a, sq(a)), frob2(a)));
 u64 T15= mul(T14, a);
 u64 acc= frob4(T15);
 for(int i= 14; i--;) acc= frob4(mul(acc, T15));
 return mul(acc, T14);
}
}  // namespace gf2_64_pclmul_window
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_window::inv;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
