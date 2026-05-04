#pragma once
// PCLMUL 直接の基本実装。 mul の比較対象なので自己完結 (= _shared/pclmul_core.hpp
// に依存しない、 _shared が将来更新されても本ファイルの mul は frozen)。
//
// 実装: P(x) = x^64 + x^4 + x^3 + x + 1 で削減する素直な PCLMUL。
//   1. clmul で 64×64 → 128 ビット多項式積
//   2. hi の各ビット (x^{64+i}) を x^{i+4} + x^{i+3} + x^{i+1} + x^i に置換 → lo に xor
//   3. 上位 4 bit (x^124..x^127) は左シフトで lo 外に出るので、 4-bit table で吸収
#pragma GCC optimize("O3,unroll-loops")
#include " ../../_shared/_common.hpp"
namespace gf2_64_mul_pclmul_baseline {
[[gnu::target("pclmul")]] [[gnu::always_inline]] inline __m128i clmul(u64 a, u64 b) {
 __m128i av{(long long)a, 0};
 __m128i bv{(long long)b, 0};
 return _mm_clmulepi64_si128(av, bv, 0);
}
[[gnu::target("pclmul")]] [[gnu::always_inline]] inline u64 reduce(__m128i v) {
 u64 lo= (u64)v[0], hi= (u64)v[1];
 lo^= hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4);
 static constexpr std::array<u64, 16> RED= [] {
  std::array<u64, 16> r{};
  for(int q= 0; q < 16; ++q) {
   u64 o= u64(q) ^ (u64(q) >> 1) ^ (u64(q) >> 3);
   r[q]= o ^ (o << 1) ^ (o << 3) ^ (o << 4);
  }
  return r;
 }();
 return lo ^ RED[hi >> 60];
}
[[gnu::target("pclmul")]] [[gnu::always_inline]] inline u64 mul(u64 a, u64 b) { return reduce(clmul(a, b)); }
}
struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_mul_pclmul_baseline::mul;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], bs[i]);
  return ans;
 }
};
