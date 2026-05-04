#pragma once
// 素朴な sqrt: a^(2^63) を 63 回の sq として実行。
//
// 設計:
//   - mul / sq は building block なので _shared/pclmul_core.hpp を使用 (current best)
//   - sqrt 戦略 (= 63 sq) は本ファイル内で自己完結 (= _shared::sqrt は使わない)
//
// 理由:
//   gf2-64-sqrt は sqrt 戦略を比較する problem。 _shared::sqrt は将来「より速い
//   sqrt (Frobenius byte-table linear map など)」で更新される可能性があり、 そう
//   なると本 baseline の "63 sq" の意味がなくなる。 計算ロジックは algo 内に閉じる。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#define PCLMUL_RUN [[gnu::target("pclmul")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_sqrt_pclmul_pow {
using gf2_64_pclmul::sq;

// sqrt(a) = a^(2^63) を 63 連続 sq で
[[gnu::target("pclmul")]] inline u64 sqrt_via_sq(u64 a) {
 for(int i= 0; i < 63; ++i) a= sq(a);
 return a;
}
}

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sqrt_pclmul_pow::sqrt_via_sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sqrt_via_sq(as[i]);
  return ans;
 }
};
