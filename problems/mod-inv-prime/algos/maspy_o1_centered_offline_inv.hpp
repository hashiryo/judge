#pragma once
#include "_common.hpp"
// =============================================================================
// maspy_o1_centered.hpp の INV テーブル precompute を offline batch で構築。
//
// 元の漸化式: K × (1 ceil-div + 1 mul + 1 mod) ≈ 30+ms
// offline batch: K × 3 mul + 1 modpow に置換、div を排除。
// テーブルは 2K+1 サイズ・K 中心構造 (INV[K-i] = -INV[K+i]) は維持。
// =============================================================================
struct ModInv {
 using i32 = int32_t;
 using u16 = unsigned short;
 static constexpr u32 K = 1u << 21;
 static constexpr u32 FRAC_BUCKETS = 1u << 20;
 static constexpr u32 STERN_LIMIT = 2048;

 static u32 pow_mod(u64 a, u64 e, u32 p) {
  u64 r = 1;
  while (e) { if (e & 1) r = r * a % p; a = a * a % p; e >>= 1; }
  return u32(r);
 }

 static vector<u32> run(u32 p, const vector<u32>& qs) {
  // ---- INV[K+1..K+K] を offline batch で構築 ----
  vector<u32> INV(2 * K + 1, 0);
  // Phase 1: INV[K+i] に i! を置く (in-place 利用)
  INV[K + 0] = 1;  // 0! = 1
  for (u32 i = 1; i <= K; ++i) {
   INV[K + i] = u32(u64(INV[K + i - 1]) * i % p);
  }
  // Phase 2: cur = (K!)^{-1}
  u64 cur = pow_mod(INV[K + K], u64(p) - 2, p);
  // Phase 3: i 降順で INV[K+i] = cur · (i-1)!  (= i^{-1})
  for (u32 i = K; i >= 1; --i) {
   u32 inv_i = u32(cur * INV[K + i - 1] % p);
   cur = cur * i % p;
   INV[K + i] = inv_i;
  }
  INV[K + 0] = 0;
  // 負側
  for (u32 i = 1; i <= K; ++i) {
   INV[K - i] = (INV[K + i] == 0) ? 0 : (p - INV[K + i]);
  }

  // ---- FRAC table 構築 (元実装と同じ) ----
  vector<pair<u16, u16>> FRAC(FRAC_BUCKETS + 1);
  std::vector<std::tuple<u32, u32, u32, u32>> stack;
  stack.emplace_back(0, 1, 1, 1);
  while (!stack.empty()) {
   auto [a, b, c, d] = stack.back(); stack.pop_back();
   if (b + d < STERN_LIMIT) {
    stack.emplace_back(a + c, b + d, c, d);
    stack.emplace_back(a, b, a + c, b + d);
    continue;
   }
   u32 s = u32(u64(a) * p / (1024u * b));
   u32 t = u32(u64(c) * p / (1024u * d));
   if (s <= FRAC_BUCKETS) FRAC[s] = {u16(a), u16(b)};
   if (t <= FRAC_BUCKETS) FRAC[t] = {u16(c), u16(d)};
   const u32 a_min = std::min(a, c), b_min = std::min(b, d);
   for (u32 i = s + 1; i < t && i <= FRAC_BUCKETS; ++i) {
    FRAC[i] = {u16(a_min), u16(b_min)};
   }
  }

  // ---- per-query ----
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 x : qs) {
   if (x == 0) { ans.push_back(u32(-1)); continue; }
   if (x == 1) { ans.push_back(1); continue; }
   const auto [a, b] = FRAC[x >> 10];
   const u32 t = x * b - u32(a) * p;
   const u32 idx = K + t;
   ans.push_back(u32(u64(INV[idx]) * b % p));
  }
  return ans;
 }
};
