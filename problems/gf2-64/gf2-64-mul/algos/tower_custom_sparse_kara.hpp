#pragma once
// 自前塔 (sparse 版) + Karatsuba (9 mul) で乗算。
// schoolbook 版の sparse (= tower_custom_sparse.hpp) との比較対象。
#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/tower_custom_sparse.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  gf2_64_tower_custom_sparse::init_f16_tables();
  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   ans[i] = gf2_64_tower_custom_sparse::mul_via_tower_karatsuba(as[i], bs[i]);
  }
  return ans;
 }
};
