#pragma once
// 自前 (= 非 Nimber) 塔表現 + Karatsuba (9 mul) で乗算を高速化した版。
// schoolbook 版 (tower_custom.hpp) の 16 → 9 mul で ~1/2 の演算量。
#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/tower_custom.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  gf2_64_tower_custom::init_f16_tables();

  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   ans[i] = gf2_64_tower_custom::mul_via_tower_karatsuba(as[i], bs[i]);
  }
  return ans;
 }
};
