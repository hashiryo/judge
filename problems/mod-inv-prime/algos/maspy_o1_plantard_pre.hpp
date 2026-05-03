#pragma once
#include "_common.hpp"
// =============================================================================
// maspy_o1_plantard.hpp の per-query を plantard_pre.hpp 流の preconditioned
// 乗算で詰めるバリアント。
//
// 元の per-query: M.reduce(u64(den) * inv_M[anum])  ≈ 3 mul
// preconditioned: M.mul_pre(den, inv_M_pre[anum])   ≈ 2 mul
//   ここで inv_M_pre[anum] = inv_M[anum] * iv mod 2^64 (u64)
//
// 期待: per-query 1 mul 削減 × 1M クエリ ≈ 0.3ms 短縮
// 代償: テーブル 5MB (u32×1.3M) → 10MB (u64×1.3M)、cache miss 増加リスク
// precompute: raw 漸化式に加え inv_M_pre 構築ループ (1.3M × 1 mul) を追加
//
// === 実測結果 (T=1M, x86_64): negative result ===
// per-query: 5.70ms → 9.14ms (+60%!)
//   anum がほぼランダム (Farey lookup 結果) で prefetcher が効かず、
//   要素サイズ倍増分そのまま LL cache miss に化けた。1 mul 削減効果は埋もれた。
// precompute: 6.47ms → 7.70ms (+1.2ms) — to_pre 追加分。
// → preconditioned lookup は「テーブルが小さい (L2 内に収まる) かつ
//   ランダムアクセス」のケース以外では基本損。教訓的サンプルとして残す。
// =============================================================================
struct ModInv {
 using i32 = int32_t;
 static constexpr u32 MAGIC1 = 1000000;
 static constexpr u32 MAGIC2 = 1300000;

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
  // preconditioned: r_pre = r * iv mod 2^64
  inline u64 to_pre(u32 r) const { return u64(r) * iv; }
  inline u32 mul_pre(u32 l, u64 r_pre) const {
   const u64 t = u64(l) * r_pre;
   return u32((u128(t | u32(-1)) * mod) >> 64);
  }
 };

 static vector<u32> run(u32 p, const vector<u32>& qs) {
  MP M(p);
  // 1) inv_lookup を raw 形で構築 (line recurrence)
  vector<u32> inv_M(MAGIC2 + 1, 0);
  inv_M[1] = 1;
  for (u32 i = 2; i <= MAGIC2; ++i) {
   inv_M[i] = u32(p - u64(p / i) * inv_M[p % i] % p);
   if (inv_M[i] == p) inv_M[i] = 0;
  }
  // 2) raw → preconditioned (Montgomery·iv = raw·R·iv mod 2^64) を u64 で保持
  // mul_pre(den, x_pre) = den · x mod p (raw 出力) になるためには
  //   x_pre = (Montgomery 形 x) · iv = (x·R) · iv mod 2^64 が必要。
  // つまり inv_M_pre[i] = M.set(inv_M_raw[i]) を Montgomery 形に変換した後 to_pre。
  vector<u64> inv_M_pre(MAGIC2 + 1, 0);
  for (u32 i = 1; i <= MAGIC2; ++i) {
   inv_M_pre[i] = M.to_pre(M.set(inv_M[i]));
  }

  // 3) Farey table
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

  // 4) per-query: 1 mul_pre (= 2 mul) で raw 値が出る
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
   u32 r = M.mul_pre(den, inv_M_pre[anum]);
   if (num < 0) r = p - r;
   if (r >= p) r -= p;
   ans.push_back(r);
  }
  return ans;
 }
};
