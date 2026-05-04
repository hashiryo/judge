#pragma once
// 素朴な binary exponentiation: a^e を 二進展開で sq + (条件付き mul)。
//
// 設計:
//   - mul / sq は building block なので _shared/pclmul_core.hpp を使用 (current best)
//   - pow 戦略 (= 単純 binary) は本ファイル内で自己完結 (= _shared::pow は使わない)
//
// 理由:
//   gf2-64-pow は pow 戦略を比較する problem。 _shared::pow が将来「より速い pow
//   戦略 (window 法 / subfield split / Frobenius など)」で更新された場合、 本
//   baseline の "naive binary" の意味がなくなる。 pow 計算ロジックは algo 内に閉じる。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_pclmul_binary {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
inline u64 pow_binary(u64 a, u64 e) {
 u64 res= 1;
 while(e) {
  if(e & 1) res= mul(res, a);
  a= sq(a);
  e>>= 1;
 }
 return res;
}
}
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pclmul_binary::pow_binary;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow_binary(as[i], es[i]);
  return ans;
 }
};
