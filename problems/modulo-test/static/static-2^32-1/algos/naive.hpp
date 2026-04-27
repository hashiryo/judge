#pragma once
#include "_common.hpp"
struct MP {  // mod < 2^32
 u32 mod;
 constexpr MP(): mod(0) {}
 constexpr MP(u32 m): mod(m) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return u64(l) * r % mod; }
 constexpr inline u32 set(u32 n) const { return n; }
 constexpr inline u32 get(u32 n) const { return n; }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const {
  l+= r;
  if(l < r) l-= mod;
  if(l >= mod) l-= mod;
  return l;
 }
 constexpr inline u32 diff(u32 l, u32 r) const { return plus(l, mod - r); }
};
