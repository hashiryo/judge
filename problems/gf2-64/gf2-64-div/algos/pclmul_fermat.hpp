#pragma once
// Fermat inverse: b^{-1} = b^(2^64 - 2)、 二進展開 pow で 63 muls + 64 sqs。
//
// 設計:
//   - mul / sq は building block なので _shared/pclmul_core.hpp を使用 (current best)
//   - inv 戦略 (= Fermat) は本ファイル内で自己完結 (= _shared::inv は使わない)
//
// 理由:
//   gf2-64-div は inv 戦略を比較する problem。 _shared::inv が将来「より速い inv
//   戦略 (Itoh-Tsujii や norm-based) で更新された場合、 本 baseline の "Fermat
//   strategy" の意味がなくなってしまう。 そのため inv 計算ロジックは algo 内に閉じる。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/pclmul_core.hpp"

namespace gf2_64_div_pclmul_fermat {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// Fermat inv: a^(2^64 - 2)。 二進展開、 各 bit で sq + (条件付き mul)。
[[gnu::target("pclmul")]] inline u64 inv_fermat(u64 a) {
 u64 res= 1;
 u64 e= ~u64(1);  // = 2^64 - 2
 while(e) {
  if(e & 1) res= mul(res, a);
  a= sq(a);
  e>>= 1;
 }
 return res;
}
}

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_div_pclmul_fermat::inv_fermat;
  using gf2_64_pclmul::mul;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv_fermat(bs[i]));
  return ans;
 }
};
