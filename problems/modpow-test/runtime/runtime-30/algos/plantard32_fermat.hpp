#pragma once
#include "_common.hpp"
// 素数 mod 限定の最適化: フェルマーの小定理で a^e ≡ a^(e mod (p-1)) (a ≠ 0) を使い、
// pow ループに入る前に exponent を mod (p-1) で削減する変種。
//
// b の bit 数が log2(p-1) より大きい (例: bbits=32, p≈2^30) と短縮効果あり。
// b ≤ p-1 のケースでは無効 (むしろ % のオーバーヘッドだけ損)。
//
// 注: 本ファイルは「mod が素数」を仮定。非素数 mod では誤った結果を返す。
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
 inline u32 pow(u32 base, u32 e) const {
  // フェルマー: a ≠ 0 のとき a^(p-1) ≡ 1。e %= p-1 で短縮可。
  // base = 0 のときは reduce 不要 (e=0 → 1, e>0 → 0)。
  if (e >= mod_minus_1) e %= mod_minus_1;
  u32 r = set(1);
  while (e) {
   if (e & 1) r = mul(r, base);
   base = mul(base, base);
   e >>= 1;
  }
  return r;
 }
};
