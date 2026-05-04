#pragma once
// sq を 4 回繰り返して a^16 を計算 (= current best sq の 4 回適用)。
//
// _shared/sq.hpp の current best sq を使う (= mul/sq は building block)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/sq.hpp"

namespace gf2_64_frob16_sq_chain {
using gf2_64_pclmul::sq;
inline u64 frob16(u64 a) {
 for(int i= 0; i < 4; ++i) a= sq(a);
 return a;
}
}

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_frob16_sq_chain::frob16;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= frob16(as[i]);
  return ans;
 }
};
