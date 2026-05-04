#pragma once
// PCLMUL + PDEP 二乗 (windowなし、binary exp のまま sq だけ高速化)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif
namespace gf2_64_pow_pdep_simple {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
[[gnu::target("pclmul")]] u64 pow(u64 a, u64 e) {
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
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_pdep_simple::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
