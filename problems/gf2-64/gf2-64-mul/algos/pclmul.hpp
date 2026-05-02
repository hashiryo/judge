#pragma once
// PCLMUL を直接使う基本実装。 _shared/pclmul_core.hpp を共有。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#define PCLMUL_RUN [[gnu::target("pclmul")]]
#else
#define PCLMUL_RUN
#endif

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = mul(as[i], bs[i]);
  return ans;
 }
};
