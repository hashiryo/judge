#pragma once
// =============================================================================
// GF(2^64) ≅ GF(2)[x] / (x^64 + x^4 + x^3 + x + 1) の基本演算 (PCLMUL 経路)。
// 5 つの実験 problem (gf2-64-mul/div/pow/sqrt/log) で共有する low-level 実装。
// 各 problem の algos/pclmul.hpp 系から #include "../../_shared/pclmul_core.hpp"
// で取り込む想定。
// =============================================================================
#include <array>
#include <cstdint>

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#elif defined(USE_SIMDE)
#include <simde/x86/sse2.h>
#include <simde/x86/clmul.h>
#else
#error "pclmul_core.hpp: requires PCLMUL (x86 native or SIMDe)."
#endif

namespace gf2_64_pclmul {
using u64 = unsigned long long;

// Carryless multiply 64×64 → 128 を 1 命令で。
[[gnu::always_inline]] inline __m128i clmul(u64 a, u64 b) {
 __m128i av{(long long) a, 0};
 __m128i bv{(long long) b, 0};
 return _mm_clmulepi64_si128(av, bv, 0);
}

// (lo, hi) を P(x) = x^64 + x^4 + x^3 + x + 1 で reduce → 64-bit
[[gnu::always_inline]] inline u64 reduce(__m128i v) {
 u64 lo = (u64) v[0], hi = (u64) v[1];
 // hi の bit i は x^(64+i)。x^64 ≡ x^4 + x^3 + x + 1 を最初に xor。
 lo ^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
 // 上位 4 bit (x^124..x^127) は左シフトで lo 外に出るので、
 // 事前計算した「上位 4 bit 由来の reduce 残渣」テーブルで吸収。
 static constexpr std::array<u64, 16> RED = [] {
  std::array<u64, 16> r{};
  for (int q = 0; q < 16; ++q) {
   u64 o = q ^ (q >> 1) ^ (q >> 3);
   r[q] = o ^ (o << 1) ^ (o << 3) ^ (o << 4);
  }
  return r;
 }();
 return lo ^ RED[hi >> 60];
}

[[gnu::always_inline]] inline u64 mul(u64 a, u64 b) {
 return reduce(clmul(a, b));
}

// 二乗。今は mul(a, a) と同じ展開をするだけ。
// 実験 algo (例: pdep_square.hpp) では PDEP で書き換え可能。
[[gnu::always_inline]] inline u64 sq(u64 a) {
 return mul(a, a);
}

// 累乗 (a^e) 二進展開。
inline u64 pow(u64 a, u64 e) {
 u64 res = 1;
 while (e) {
  if (e & 1) res = mul(res, a);
  a = sq(a);
  e >>= 1;
 }
 return res;
}

// Fermat-style inverse: a^(2^64 - 2)。
// より高速な Itoh-Tsujii は別 algo で試す。
inline u64 inv(u64 a) {
 // exponent = 2^64 - 2 = 0xFFFFFFFFFFFFFFFE
 return pow(a, ~(u64) 1);
}

// sqrt: Frobenius (^2) の逆 = ^(2^63)。
// 専用の linear-map 実装は別 algo で試す。
inline u64 sqrt(u64 a) {
 return pow(a, (u64) 1 << 63);
}

} // namespace gf2_64_pclmul
