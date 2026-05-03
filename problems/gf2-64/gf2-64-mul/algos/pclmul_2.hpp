#pragma once
// PCLMUL 直接の基本実装。 mul の比較対象なので自己完結 (= _shared/pclmul_core.hpp
// に依存しない、 _shared が将来更新されても本ファイルの mul は frozen)。
//
// 実装: P(x) = x^64 + x^4 + x^3 + x + 1 で削減する素直な PCLMUL。
//   1. clmul で 64×64 → 128 ビット多項式積
//   2. hi の各ビット (x^{64+i}) を x^{i+4} + x^{i+3} + x^{i+1} + x^i に置換 → lo に xor
//   3. 上位 4 bit (x^124..x^127) は左シフトで lo 外に出るので、 4-bit table で吸収
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
#error "pclmul.hpp: requires PCLMUL (x86 native or SIMDe)."
#endif
namespace gf2_64_mul_pclmul_baseline {
inline u64 mul(u64 a, u64 b) {
 static constexpr std::array<u64, 8> RED= {0, 27, 45, 54, 90, 65, 119, 108};
 __m128i v= _mm_clmulepi64_si128(_mm_cvtsi64_si128(a), _mm_cvtsi64_si128(b), 0);
 u64 h= (u64)v[1], d= h ^ (h << 1);
 return (u64)v[0] ^ RED[h >> 60] ^ d ^ (d << 3);
}
}
struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_mul_pclmul_baseline::mul;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], bs[i]);
  return ans;
 }
};
