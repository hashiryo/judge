#pragma once
#include "_common.hpp"
// =============================================================================
// Maspy 氏自身による公開実装 (modfast.hpp) を参考にした O(1) inv バリアント。
//
// 移植元: https://maspypy.github.io/library/mod/modfast.hpp
//
// purplesyringa 提出 (=我々の maspy_o1.hpp の元) と比べた特徴:
//   1. INV テーブルを 2K+1 サイズで K 中心: INV[K+i] = i^{-1} (i ∈ [-K, K])
//      → per-query で abs() と (num < 0) 分岐がゼロになる
//   2. Stern-Brocot 木 BFS で FRAC を構築 (Farey 再帰より単純)
//
// per-query (~5 cycle、ours の ~8 cycle より短い):
//   auto [a, b] = FRAC[x >> 10];
//   u32 t = x * b - a * p;          // u32 wrap で「負」も index 化
//   return INV[K + t] * u64(b) % p; // 分岐なし
//
// パラメータ (maspy 氏 modfast.hpp と同じ):
//   K = 2^21 = 2097152                  (INV size = 2K+1 = ~4M, ~16 MB)
//   FRAC bucket = x >> 10, 約 1M entries (= p/1024) サイズ
//   Stern-Brocot 上限 b+d < 2048
//
// 注: p = 998244353 < 2^30 を仮定。
// =============================================================================
struct ModInv {
 using i32 = int32_t;
 using u16 = unsigned short;
 static constexpr u32 K = 1u << 21;          // = 2097152
 static constexpr u32 FRAC_BUCKETS = 1u << 20;  // = 1048576 (≈ p/1024)
 static constexpr u32 STERN_LIMIT = 2048;    // b+d がこれ以上で葉

 static vector<u32> run(u32 p, const vector<u32>& qs) {
  // ---- INV table 構築 (Maspy 流の漸化式) ----
  // INV[K + i] = i^{-1} mod p for i ∈ [1, K]
  // 漸化式: q = ceil(p/i)、INV[K+i] = INV[K + i*q - p] · q mod p
  //   検証: i · INV[K+i] = i · q · INV[K + iq - p] = (iq - p + p) · ...
  //                    ≡ (iq - p) · INV[iq-p] ≡ 1 (mod p)
  vector<u32> INV(2 * K + 1, 0);
  INV[K + 1] = 1;
  for (u32 i = 2; i <= K; ++i) {
   u64 q = (u64(p) + i - 1) / i;  // ceil(p/i)
   u32 iq_minus_p = u32(u64(i) * q - p);
   INV[K + i] = u32(u64(INV[K + iq_minus_p]) * q % p);
  }
  // 負側: INV[K - i] = p - INV[K + i]
  for (u32 i = 1; i <= K; ++i) {
   INV[K - i] = (INV[K + i] == 0) ? 0 : (p - INV[K + i]);
  }

  // ---- FRAC table 構築 (Stern-Brocot 木 BFS) ----
  // FRAC[bucket] = (a, b) where a/b ≈ x/p for x in this bucket
  vector<pair<u16, u16>> FRAC(FRAC_BUCKETS + 1);
  // BFS スタック: (a, b, c, d) で a/b と c/d の中項を再帰的に挿入
  // 初期: 0/1 と 1/1 (Stern-Brocot 木のルート 2 ペア)
  std::vector<std::tuple<u32, u32, u32, u32>> stack;
  stack.emplace_back(0, 1, 1, 1);
  while (!stack.empty()) {
   auto [a, b, c, d] = stack.back(); stack.pop_back();
   if (b + d < STERN_LIMIT) {
    // 中項 (a+c)/(b+d) を挿入、左右に分割再帰
    stack.emplace_back(a + c, b + d, c, d);
    stack.emplace_back(a, b, a + c, b + d);
    continue;
   }
   // 葉に到達 — bucket index に書き込み
   u32 s = u32(u64(a) * p / (1024u * b));
   u32 t = u32(u64(c) * p / (1024u * d));
   if (s <= FRAC_BUCKETS) FRAC[s] = {u16(a), u16(b)};
   if (t <= FRAC_BUCKETS) FRAC[t] = {u16(c), u16(d)};
   const u32 a_min = std::min(a, c), b_min = std::min(b, d);
   for (u32 i = s + 1; i < t && i <= FRAC_BUCKETS; ++i) {
    FRAC[i] = {u16(a_min), u16(b_min)};
   }
  }

  // ---- Per-query: O(1) inv ----
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 x : qs) {
   if (x == 0) { ans.push_back(u32(-1)); continue; }
   if (x == 1) { ans.push_back(1); continue; }
   const auto [a, b] = FRAC[x >> 10];
   const u32 t = x * b - u32(a) * p;  // u32 wrap で signed 値を代用
   // K + t を u32 で計算: t が「負」(= u32 wrap で大きい値) でも、
   //   K (= 2^21) を加えると 2^21 + (2^32 + signed_t) = 2^32 + (2^21 + signed_t)
   //   u32 で取ると K + signed_t (signed_t が負なら K + signed_t < K)
   // |signed_t| ≤ K となるよう Stern-Brocot 上限を設計しているので index 範囲内
   const u32 idx = K + t;  // u32 wrap で正しい index
   ans.push_back(u32(u64(INV[idx]) * b % p));
  }
  return ans;
 }
};
