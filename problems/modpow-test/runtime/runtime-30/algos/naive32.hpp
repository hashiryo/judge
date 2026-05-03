#pragma once
#include "_common.hpp"
// 標準的な u32 mod 演算 + binary pow。素朴ベースライン。
struct MP {  // mod < 2^32
 u32 mod;
 constexpr MP(): mod(0) {}
 constexpr MP(u32 m): mod(m) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return u32(u64(l) * r % mod); }
 constexpr inline u32 set(u32 n) const { return n; }
 constexpr inline u32 get(u32 n) const { return n; }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l += r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l -= r, l >> 31 ? l + mod : l; }
 inline u32 pow(u32 base, u32 e) const {
  u32 r = set(1);
  while (e) {
   if (e & 1) r = mul(r, base);
   base = mul(base, base);
   e >>= 1;
  }
  return r;
 }
};
