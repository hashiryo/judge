#pragma once
#include "_common.hpp"
// Plantard 32-bit + base 側を preconditioned 形 (base * iv mod 2^64) で保持し、
// pow ループの r *= base 部分で 1 mul/call 削る変種。
// 詳細は mod-inv-prime/algos/plantard_pre.hpp と同じ手筋。
struct MP {  // mod < 2^32/phi
 u32 mod;
 u32 r2;
 u64 iv;
 constexpr MP(): mod(0), r2(0), iv(0) {}
 constexpr MP(u32 m): mod(m), r2(u32(-u128(m) % m)), iv(inv_(m)) {}
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
