#pragma once
#include "_common.hpp"
// https://gmplib.org/~tege/division-paper.pdf
struct MP {  // 2^31 <= mod < 2^32
 u64 mod;
 constexpr MP(): mod(0), v(0) {}
 constexpr MP(u64 m): mod(m), v(u128(-1) / m) {}
 constexpr inline u64 mul(u64 l, u64 r) const { return rem(u128(l) * r); }
 constexpr inline u64 set(u64 n) const { return n; }
 constexpr inline u64 get(u64 n) const { return n; }
 constexpr inline u64 norm(u64 n) const { return n; }
 constexpr inline u64 plus(u64 l, u64 r) const {
  l+= r;
  if(l < r) l-= mod;
  if(l >= mod) l-= mod;
  return l;
 }
 constexpr inline u64 diff(u64 l, u64 r) const { return plus(l, mod - r); }
private:
 u64 v;
 constexpr inline u64 rem(u128 n) const {
  u128 q= u128(v) * u64(n >> 64) + n;
  u64 q1= u64(q >> 64) + 1;
  u64 r= n - q1 * mod;
  if(r > u64(q)) r+= mod;
  if(r >= mod) r-= mod;
  return r;
 }
};
