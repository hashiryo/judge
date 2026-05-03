#pragma once
// sq(a) = mul(a, a) を PCLMUL で計算 (素朴ベースライン)。
// 自己完結 (= _shared/pclmul_core.hpp に依存しない、 sq の比較対象なので)。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#endif
#include "_common.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul")]]
#elif defined(USE_SIMDE)
#include <simde/x86/sse2.h>
#include <simde/x86/clmul.h>
#define PCLMUL_RUN
#else
#error "pclmul_mul.hpp: requires PCLMUL (x86 native or SIMDe)."
#endif

namespace gf2_64_sq_pclmul_mul {

[[gnu::always_inline]] inline u64 sq(u64 a) {
 __m128i av{(long long) a, 0};
 __m128i v = _mm_clmulepi64_si128(av, av, 0);
 u64 lo = (u64) v[0], hi = (u64) v[1];
 lo ^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
 static constexpr std::array<u64, 16> RED = [] {
  std::array<u64, 16> r{};
  for (int q = 0; q < 16; ++q) {
   u64 o = u64(q) ^ (u64(q) >> 1) ^ (u64(q) >> 3);
   r[q] = o ^ (o << 1) ^ (o << 3) ^ (o << 4);
  }
  return r;
 }();
 return lo ^ RED[hi >> 60];
}

}

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_pclmul_mul::sq;
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = sq(as[i]);
  return ans;
 }
};
