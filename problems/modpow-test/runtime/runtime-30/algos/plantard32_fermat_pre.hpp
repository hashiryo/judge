#pragma once
#include "_common.hpp"
// Fermat 削減 (e %= p-1) + preconditioned base の合わせ技。
//   - 高 bbits 領域で iteration 削減
//   - 残ったループでも r *= base を mul_pre で 1 mul 削減
// 期待: bbits >= log2(p-1) + α で plantard32_fermat / plantard32_pre 双方を上回る
//
// === 実測結果 (n=5M, p=998244353): negative result ===
// heavy_00 (bbits=30): 548ms vs plantard32_pre 506ms (+8%、両方の利点を失う)
// heavy_01 (bbits=32): 557ms vs plantard32_fermat 531ms (+5%)
// 推察: _pre 版のループは `if (!(e>>=1)) break;` 構造で分岐予測 / scheduler に
//   敏感。divq (e %= mod_minus_1) を頭に追加するとクリティカルパスが伸び、
//   ループ最適化が崩れる。「per-iter mul 削減」と「iteration 数削減」を独立に
//   合算できるはず、という素朴な仮定が成り立たない例。
struct MP {  // mod < 2^32/phi, mod は素数
 u32 mod;
 u32 mod_minus_1;
 u32 r2;
 u64 iv;
 constexpr MP(): mod(0), mod_minus_1(0), r2(0), iv(0) {}
 constexpr MP(u32 m): mod(m), mod_minus_1(m - 1), r2(u32(-u128(m) % m)), iv(inv_(m)) {}
 static constexpr u64 inv_(u64 n, int e = 6, u64 x = 1) {
  return e ? inv_(n, e - 1, x * (2 - x * n)) : x;
 }
 constexpr inline u32 reduce(u64 w) const { return u32((u128((w * iv) | u32(-1)) * mod) >> 64); }
 constexpr inline u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return mul(n, r2); }
 constexpr inline u32 get(u32 n) const { u32 v = reduce(n); return v >= mod ? v - mod : v; }
 constexpr inline u32 norm(u32 n) const { return n >= mod ? n - mod : n; }
 inline u64 to_pre(u32 r) const { return u64(r) * iv; }
 inline u32 mul_pre(u32 l, u64 r_pre) const {
  const u64 t = u64(l) * r_pre;
  return u32((u128(t | u32(-1)) * mod) >> 64);
 }
 inline u32 pow(u32 base, u32 e) const {
  if (e >= mod_minus_1) e %= mod_minus_1;
  u32 r = set(1);
  if (!e) return r;
  u64 base_pre = to_pre(base);
  for (;;) {
   if (e & 1) r = mul_pre(r, base_pre);
   if (!(e >>= 1)) break;
   base = mul_pre(base, base_pre);
   base_pre = to_pre(base);
  }
  return r;
 }
};
