#pragma once
#include "_common.hpp"
// =============================================================================
// Maspy 流 O(1) inv:
//   x^{-1} を x ≡ a/b mod p (|a| ≤ 2p/n, b ≤ n) の rational 分解で求める。
//   x^{-1} = b · a^{-1} (or -b · |a|^{-1} when a < 0)
//   |a| ≤ 2p/n の範囲の inv を事前計算 (linear recurrence で O(1) per entry)。
//
// 参考: https://maspypy.com/o1-mod-inv-mod-pow
//
// パラメータ (n = p^{1/3} 周辺で最適、purplesyringa の log 提出と同じ):
//   MAGIC1 = 1000000   // Farey buckets ([0, 1) を 10^6 分割)
//   MAGIC2 = 1300000   // |a| ≤ 2p/n の上限 ≈
//
// per-query: 1 bucket lookup + 1 inv lookup + 1 mul (+条件付き sub) ≈ 10 cycle
// precompute: ~50 ms (linear recurrence for inv table + Farey 構築)
// =============================================================================
struct ModInv {
 using i32 = int32_t;
 static constexpr u32 MAGIC1 = 1000000;
 static constexpr u32 MAGIC2 = 1300000;

 static vector<u32> run(u32 p, const vector<u32>& qs) {
  // ---- 事前計算 ----
  // inv_lookup[i] = i^{-1} mod p, for i ∈ [1, MAGIC2]
  // 線型漸化式: inv[i] = -(p/i) · inv[p%i] mod p
  vector<u32> inv_lookup(MAGIC2 + 1, 0);
  inv_lookup[1] = 1;
  for (u32 i = 2; i <= MAGIC2; ++i) {
   u32 q = p / i, r = p % i;
   // r < i なので inv_lookup[r] は既に計算済 (r < MAGIC2 自明)
   inv_lookup[i] = u32(p - u64(q) * inv_lookup[r] % p);
   if (inv_lookup[i] == p) inv_lookup[i] = 0;
  }

  // Farey table: bucket → (numerator, denominator) packed in u32
  // 上位 16 bit: c (numerator), 下位 16 bit: b (denominator)
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

  // ---- per-query ----
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
   // x = a/b で a = num, b = den ⇒ x^{-1} = b · a^{-1}
   u32 r = u32(u64(den) * inv_lookup[anum] % p);
   if (num < 0) r = (r == 0) ? 0 : p - r;
   ans.push_back(r);
  }
  return ans;
 }
};
