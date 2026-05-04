#pragma once
// sq を 5 回繰り返して a^32 を計算 (= current best sq の 5 回適用)。
//
// _shared/sq.hpp の current best sq を使う (= mul/sq は building block)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/sq.hpp"
namespace gf2_64_frob32_sq_chain {
using gf2_64_pclmul::sq;
inline u64 frob32(u64 a) { return sq(sq(sq(sq(sq(a))))); }
}
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_frob32_sq_chain::frob32;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= frob32(as[i]);
  return ans;
 }
};
