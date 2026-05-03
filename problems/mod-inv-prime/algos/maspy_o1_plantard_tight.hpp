#pragma once
#include "_common.hpp"
// =============================================================================
// maspy_o1_plantard.hpp の MAGIC2 を理論下限 (= p/√MAGIC1 ≈ 998244) に揃えた版。
//
// 元実装 purplesyringa の MAGIC2 = 1.3M は安全マージン 30% 込みだが、実測上
// 不要 (random_large_00 全 1M クエリで一致) で、inv table を 5.2MB → 4MB に
// 縮めた分 cache footprint が改善し速くなる。
//
// 採用: MAGIC1 = 1M, MAGIC2 = 1M (理論下限ぴったり)
//   √MAGIC1 = 1000, MAGIC2 × √MAGIC1 = 10^9 ≳ p (998244353) ✓
//   memory: Farey 4MB + inv 4MB = 8MB (元 9.2MB)
// =============================================================================
struct ModInv {
 using i32 = int32_t;
 static constexpr u32 MAGIC1 = 1000000;
 static constexpr u32 MAGIC2 = 1000000;

 struct MP {
  u32 mod;
  u32 r2;
  u64 iv;
  constexpr MP(u32 m): mod(m), r2(u32(-u128(m) % m)), iv(inv_(m)) {}
  static constexpr u64 inv_(u64 n, int e = 6, u64 x = 1) {
   return e ? inv_(n, e - 1, x * (2 - x * n)) : x;
  }
  constexpr u32 reduce(u64 w) const { return u32((u128((w * iv) | u32(-1)) * mod) >> 64); }
  constexpr u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
  constexpr u32 set(u32 n) const { return mul(n, r2); }
  constexpr u32 get(u32 n) const { return reduce(n); }
 };

 static vector<u32> run(u32 p, const vector<u32>& qs) {
  MP M(p);
  vector<u32> inv_M(MAGIC2 + 1, 0);
  inv_M[1] = 1;
  for (u32 i = 2; i <= MAGIC2; ++i) {
   inv_M[i] = u32(p - u64(p / i) * inv_M[p % i] % p);
   if (inv_M[i] == p) inv_M[i] = 0;
  }
  for (u32 i = 1; i <= MAGIC2; ++i) inv_M[i] = M.set(inv_M[i]);

  vector<u32> farey_lookup(MAGIC1, 0);
  auto farey_rec = [&](auto& self, u32 f1, u32 f2, u32 x, u32 y) -> void {
   u32 f3 = f1 + f2;
   u32 lo = (((u64) p * (f3 >> 16) - MAGIC2) * MAGIC1 - 1) / ((u64) p * (f3 & 0xffff)) + 1;
   u32 hi = (((u64) p * (f3 >> 16) + MAGIC2) * MAGIC1) / ((u64) p * (f3 & 0xffff));
   lo = std::max(lo, x); hi = std::min(hi, y);
   if (x < lo) self(self, f1, f3, x, lo);
   std::fill(farey_lookup.begin() + lo, farey_lookup.begin() + hi, f3);
   if (hi < y) self(self, f3, f2, hi, y);
  };
  const u32 first_x = u64(MAGIC2) * MAGIC1 / p;
  const u32 first_y = (u64(p - MAGIC2) * MAGIC1 - 1) / (p * 2) + 1;
  std::fill(farey_lookup.begin(), farey_lookup.begin() + first_x, 1u);
  farey_rec(farey_rec, 1u, 0x10002u, first_x, first_y);
  std::fill(farey_lookup.begin() + first_y, farey_lookup.begin() + MAGIC1/2, 0x10002u);
  for (u32 i = MAGIC1/2; i < MAGIC1; ++i) {
   farey_lookup[i] = (farey_lookup[MAGIC1 - 1 - i] * 0xffff0001u ^ 0xffff0000u) + 0x10000u;
  }

  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 x : qs) {
   if (x == 0) { ans.push_back(u32(-1)); continue; }
   if (x == 1) { ans.push_back(1); continue; }
   const u32 bucket = u64(x) * MAGIC1 / p;
   const u32 frac = farey_lookup[bucket];
   const u32 den = frac & 0xffff;
   const i32 num = (i32) (den * x - (frac >> 16) * p);
   const u32 anum = (u32) std::abs(num);
   u32 r = M.reduce(u64(den) * inv_M[anum]);
   if (num < 0) r = p - r;
   if (r >= p) r -= p;
   ans.push_back(r);
  }
  return ans;
 }
};
