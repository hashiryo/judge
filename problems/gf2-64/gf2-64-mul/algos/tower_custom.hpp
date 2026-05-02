#pragma once
// 自前 (= 非 Nimber) 塔表現での乗算。
//   F_{2^16} = GF(2)[s]/(s^16+s^12+s^3+s+1) (primitive)
//   F_{2^64} = F_{2^16}[ω]/q(ω)   (q は ω の min poly)
// 入出力は poly 基底 u64 (他 algo と互換)、内部で塔基底に変換。
#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/tower_custom.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  // 起動時に F_{2^16} log/exp テーブル構築 (~1 ms)
  gf2_64_tower_custom::init_f16_tables();

  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   ans[i] = gf2_64_tower_custom::mul_via_tower(as[i], bs[i]);
  }
  return ans;
 }
};
