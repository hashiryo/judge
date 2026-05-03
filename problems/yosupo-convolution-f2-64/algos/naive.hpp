#pragma once
#include "_common.hpp"
// O(nm) 素朴: pclmul mul を nested loop で。 大きい n, m では TLE する想定。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#include <immintrin.h>
#endif

namespace conv_f2_64_naive {
// F_{2^64} = F_2[x] / (x^64 + x^4 + x^3 + x + 1)
// 削減多項式の table (上位 4 bit reduction 残渣)
inline const std::array<u64, 16> RED = [] {
 std::array<u64, 16> r{};
 for (int q = 0; q < 16; ++q) {
  u64 o = (u64) q ^ ((u64) q >> 1) ^ ((u64) q >> 3);
  r[q] = o ^ (o << 1) ^ (o << 3) ^ (o << 4);
 }
 return r;
}();

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
[[gnu::target("pclmul")]]
#endif
inline u64 mul(u64 a, u64 b) {
 __m128i av{(long long) a, 0};
 __m128i bv{(long long) b, 0};
 __m128i v = _mm_clmulepi64_si128(av, bv, 0);
 u64 lo = (u64) v[0], hi = (u64) v[1];
 lo ^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
 return lo ^ RED[hi >> 60];
}
}

struct Solver {
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
 [[gnu::target("pclmul")]]
#endif
 static std::vector<u64> run(int n, int m, const std::vector<u64>& a, const std::vector<u64>& b) {
  std::vector<u64> c(n + m - 1, 0);
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < m; ++j) {
    c[i + j] ^= conv_f2_64_naive::mul(a[i], b[j]);
   }
  }
  return c;
 }
};
