#pragma once
// 自前塔 (sparse 版): ω を探索して q をなるべく疎にした実装。schoolbook 16 mul。
//   現状の選択: ω = 0x20ed58e363634020、q = y^4 + 0x24a8 y^3 + 0xf994 y^2 + 0 y + 0x4c76
//   (4-term q、a_1 = 0)。reduce で 11 muls (元 12 muls)、約 -8% mul 削減。
#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/tower_custom_sparse.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  gf2_64_tower_custom_sparse::init_f16_tables();
  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   ans[i] = gf2_64_tower_custom_sparse::mul_via_tower(as[i], bs[i]);
  }
  return ans;
 }
};
