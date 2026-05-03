#pragma once
// PDEP-based squaring:
//   GF(2) では (sum b_i x^i)^2 = sum b_i x^{2i} (cross term ゼロ)。
//   bit i → bit 2i の spread は PDEP で 1 命令:
//     lo_poly = pdep(a & 0xFFFFFFFF, 0x5555...5555)  ← bit 0..31 → 0,2,..,62
//     hi_poly = pdep(a >> 32, 0x5555...5555)         ← bit 32..63 → 0,2,..,62 (= 64,66,..,126)
//   その後 P(x) で reduce (= 標準 reduction)。
//
// BMI2 が無い環境 (古い x86, ARM) ではビット並びを spread する loop で fallback。
//
// PCLMUL 比較:
//   - PCLMUL sq: 1 PCLMUL (~5 cyc) + reduce
//   - PDEP sq:   2 PDEP    (~3 cyc each, 並列) + reduce
//   PCLMUL はレイテンシ長め、 PDEP は並列性高い。 throughput は PDEP がやや勝つ可能性
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("bmi2")
#endif
#include "_common.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define BMI2_RUN [[gnu::target("bmi2")]]
#define HAS_PDEP 1
#elif defined(__BMI2__)
#define BMI2_RUN
#define HAS_PDEP 1
#else
#define BMI2_RUN
#define HAS_PDEP 0
#endif

namespace gf2_64_sq_pdep {

[[gnu::always_inline]] inline u64 spread_bits(u32 a) {
#if HAS_PDEP
 return _pdep_u64(u64(a), 0x5555555555555555ull);
#else
 // fallback: bit interleave for 32-bit input
 u64 x = a;
 x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x = (x | (x <<  8)) & 0x00FF00FF00FF00FFull;
 x = (x | (x <<  4)) & 0x0F0F0F0F0F0F0F0Full;
 x = (x | (x <<  2)) & 0x3333333333333333ull;
 x = (x | (x <<  1)) & 0x5555555555555555ull;
 return x;
#endif
}

BMI2_RUN inline u64 sq(u64 a) {
 const u64 lo = spread_bits(u32(a));
 const u64 hi = spread_bits(u32(a >> 32));
 // reduce P(x) = x^64 + x^4 + x^3 + x + 1
 u64 r = lo;
 r ^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
 static constexpr std::array<u64, 16> RED = [] {
  std::array<u64, 16> rt{};
  for (int q = 0; q < 16; ++q) {
   u64 o = u64(q) ^ (u64(q) >> 1) ^ (u64(q) >> 3);
   rt[q] = o ^ (o << 1) ^ (o << 3) ^ (o << 4);
  }
  return rt;
 }();
 return r ^ RED[hi >> 60];
}

}

struct GF2_64Op {
 BMI2_RUN static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_pdep::sq;
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = sq(as[i]);
  return ans;
 }
};
