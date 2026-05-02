#pragma once
// PCLMUL + binary exponentiation。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pclmul::pow;
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = pow(as[i], es[i]);
  return ans;
 }
};
