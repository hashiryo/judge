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
#include "../../_shared/_common.hpp"
namespace gf2_64_sq_pdep {
[[gnu::target("bmi2")]] [[gnu::always_inline]] inline u64 spread_bits(u32 a) {
#if HAS_PDEP
 return _pdep_u64(u64(a), 0x5555555555555555ull);
#else
 // fallback: bit interleave for 32-bit input
 u64 x= a;
 x= (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x= (x | (x << 8)) & 0x00FF00FF00FF00FFull;
 x= (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
 x= (x | (x << 2)) & 0x3333333333333333ull;
 x= (x | (x << 1)) & 0x5555555555555555ull;
 return x;
#endif
}
[[gnu::target("bmi2")]] inline u64 sq(u64 a) {
 u64 h= spread_bits(u32(a >> 32)), d= h ^ (h << 1);
 return spread_bits(u32(a)) ^ ((u8[]){0, 27, 45, 54, 90, 65, 119, 108})[h >> 60] ^ d ^ (d << 3);
}
}
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_pdep::sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sq(as[i]);
  return ans;
 }
};
