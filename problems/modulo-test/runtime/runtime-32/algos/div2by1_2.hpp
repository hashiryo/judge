#pragma once
#include "_common.hpp"
// https://gmplib.org/~tege/division-paper.pdf
struct MP {  // 2^31 <= mod < 2^32
 u32 mod;
 constexpr MP(): mod(0), v(0) {}
 constexpr MP(u32 m): mod(m), v(u64(-1) / m) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return n; }
 constexpr inline u32 get(u32 n) const { return n; }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u64 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u64 l, u32 r) const { return l-= r, l >> 63 ? l + mod : l; }
private:
 u32 v;
 constexpr inline u32 rem(u64 n) const {
  u64 q= u64(v) * u32(n >> 32) + n;
  u32 q1= u32(q >> 32) + 1;
  u32 r= n - q1 * mod;
  if(r > u32(q)) r+= mod;
  if(r >= mod) r-= mod;
  return r;
 }
};
