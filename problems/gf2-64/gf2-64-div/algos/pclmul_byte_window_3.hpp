#pragma once

#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pclmul_itoh_tsujii {
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob3;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob8;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
u64 inv(u64 a) {
 // 254
 u64 TFE= sq(a);
 for(int i= 6; i--;) TFE= sq(mul(a, TFE));
 u64 TFF= mul(TFE, a);
 u64 acc= frob8(TFF);
 for(int i= 6; i--;) acc= frob8(mul(acc, TFF));
 return mul(acc, TFE);
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
