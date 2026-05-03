#pragma once
#include "_common.hpp"
// =============================================================================
// maspy_o1.hpp の inv_lookup precompute を offline batch (prefix product + 1 inv)
// で構築するバリアント。
//
// 元の漸化式版: 1.3M × (1 div + 1 mul + 1 sub) ≈ 10–15ms
// offline batch:
//   1) fact[i] = i! mod p  (1.3M × 1 mul)
//   2) inv_fact[N] = pow(fact[N], p-2) (1 modpow)
//   3) i 降順で inv[i] = inv_fact[i] · fact[i-1]、cur *= i (1.3M × 2 mul)
// 合計 1.3M × 3 mul + 1 modpow。div 排除で理論上 2-3 倍速。
// 実装は in-place: 1 配列で fact → inv に書き換え。
// =============================================================================
struct ModInv {
 using i32 = int32_t;
 static constexpr u32 MAGIC1 = 1000000;
 static constexpr u32 MAGIC2 = 1300000;

 static u32 pow_mod(u64 a, u64 e, u32 p) {
  u64 r = 1;
  while (e) { if (e & 1) r = r * a % p; a = a * a % p; e >>= 1; }
  return u32(r);
 }

 static vector<u32> run(u32 p, const vector<u32>& qs) {
  // ---- offline batch で inv_lookup[1..MAGIC2] を構築 ----
  // Phase 1: fact[i] = i! mod p (in-place で inv_lookup を fact として使用)
  vector<u32> inv_lookup(MAGIC2 + 1);
  inv_lookup[0] = 1;
  for (u32 i = 1; i <= MAGIC2; ++i) {
   inv_lookup[i] = u32(u64(inv_lookup[i - 1]) * i % p);
  }
  // Phase 2: cur = (MAGIC2!)^{-1}
  u64 cur = pow_mod(inv_lookup[MAGIC2], u64(p) - 2, p);
  // Phase 3: i 降順で inv[i] = cur · fact[i-1]; cur *= i
  // inv_lookup[i-1] (= (i-1)!) は次以降の iteration では使われないので in-place 上書き可
  for (u32 i = MAGIC2; i >= 1; --i) {
   u32 inv_i = u32(cur * inv_lookup[i - 1] % p);
   cur = cur * i % p;  // (i!)^{-1} · i = ((i-1)!)^{-1}
   inv_lookup[i] = inv_i;
  }
  inv_lookup[0] = 0;

  // ---- Farey table 構築 (元実装と同じ) ----
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

  // ---- per-query (元実装と同じ) ----
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
   u32 r = u32(u64(den) * inv_lookup[anum] % p);
   if (num < 0) r = (r == 0) ? 0 : p - r;
   ans.push_back(r);
  }
  return ans;
 }
};
