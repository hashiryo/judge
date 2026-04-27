#pragma once
#include "_common.hpp"
struct MP {  // mod < 2^64
 u64 mod;
 constexpr MP(): mod(0) {}
 constexpr MP(u64 m): mod(m) {}
 constexpr inline u64 mul(u64 l, u64 r) const { return u128(l) * r % mod; }
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
};
