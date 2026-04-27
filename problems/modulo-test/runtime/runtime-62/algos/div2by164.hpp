#pragma once
#include "_common.hpp"
// https://gmplib.org/~tege/division-paper.pdf
struct MP {  // mod < 2^32
 u64 mod;
 constexpr MP(): mod(0), s(0), d(0), v(0) {}
 constexpr MP(u64 m): mod(m), s(__builtin_clzll(m)), d(m << s), v(u128(-1) / d) {}
 constexpr inline u64 mul(u64 l, u64 r) const { return rem(u128(l) * r); }
 constexpr inline u64 set(u64 n) const { return n; }
 constexpr inline u64 get(u64 n) const { return n; }
 constexpr inline u64 norm(u64 n) const { return n; }
 constexpr inline u64 plus(u64 l, u64 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u64 diff(u64 l, u64 r) const { return l-= r, l >> 63 ? l + mod : l; }
private:
 u8 s;
 u64 d;
 u64 v;
 constexpr inline u64 rem(u128 n) const {
  n<<= s;
  u128 q= u128(v) * u64(n >> 64) + n;
  u64 q1= u64(q >> 64) + 1;
  u64 r= n - q1 * d;
  if(r > u64(q)) r+= d;
  if(r >= d) r-= d;
  return r >> s;
 }
};
