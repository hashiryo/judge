#pragma once
// 素朴 reference (golden): bit-by-bit clmul + 長除法 reduce。PCLMUL 使わない。
#include "_common.hpp"

namespace gf2_64_ref {
inline std::pair<u64, u64> clmul_loop(u64 a, u64 b) {
 u64 lo = 0, hi = 0;
 for (int i = 0; i < 64; ++i) {
  if ((b >> i) & 1) {
   lo ^= a << i;
   if (i) hi ^= a >> (64 - i);
  }
 }
 return {lo, hi};
}
inline u64 reduce_naive(u64 lo, u64 hi) {
 for (int i = 63; i >= 0; --i) {
  if ((hi >> i) & 1) {
   hi ^= u64(1) << i;
   lo ^= IRRED_LOW << i;
   if (i > 0) hi ^= IRRED_LOW >> (64 - i);
  }
 }
 return lo;
}
inline u64 mul(u64 a, u64 b) {
 auto [lo, hi] = clmul_loop(a, b);
 return reduce_naive(lo, hi);
}
} // namespace gf2_64_ref

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = gf2_64_ref::mul(as[i], bs[i]);
  return ans;
 }
};
