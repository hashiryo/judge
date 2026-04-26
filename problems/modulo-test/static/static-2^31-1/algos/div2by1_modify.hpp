#pragma once
#include "_common.hpp"
// https://gmplib.org/~tege/division-paper.pdf
struct MP {  // mod < 2^31
 u32 mod;
 constexpr MP(): mod(0), s(0), d(0), v(0) {}
 constexpr MP(u32 m): mod(m), s(__builtin_clz(m)), d(m << s), v(u64(-1) / d) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return n; }
 constexpr inline u32 get(u32 n) const { return n; }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + mod : l; }
private:
 u8 s;
 u32 d;
 u32 v;
 constexpr inline u32 rem(u64 n) const {
  n<<= s;
  u64 q= u64(v) * u32(n >> 32) + n;
  u32 q1= u32(q >> 32) + 1;
  u32 r= n - q1 * d;
  r+= d & -u32(r > u32(q));
  if(__builtin_expect(r >= d, 0)) r-= d;
  return r >> s;
 }
};
