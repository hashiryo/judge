#pragma once
#include "_common.hpp"
// https://zenn.dev/yatyou/articles/f8a3969c6e9f7b
struct MP {  // mod < 2^64/phi
 u64 mod;
 constexpr MP(): mod(0), r2(0), iv(0) {}
 constexpr MP(u64 m): mod(m), r2(R2(m)), iv(inv(m)) {}
 constexpr inline u64 mul(u64 l, u64 r) const { return reduce(u128(l) * r); }
 constexpr inline u64 set(u64 n) const { return mul(n, r2); }
 constexpr inline u64 get(u64 n) const { return n= reduce(n), n >= mod ? n - mod : n; }
 constexpr inline u64 norm(u64 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u64 plus(u64 l, u64 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u64 diff(u64 l, u64 r) const { return l-= r, l >> 63 ? l + mod : l; }
 inline u64 pow(u64 base, u64 e) const {
  u64 r = set(1);
  while (e) {
   if (e & 1) r = mul(r, base);
   base = mul(base, base);
   e >>= 1;
  }
  return r;
 }
private:
 u64 r2;
 u128 iv;
 static constexpr u128 inv(u128 n, int e= 7, u128 x= 1) { return e ? inv(n, e - 1, x * (2 - x * n)) : x; }
 constexpr inline u64 reduce(const u128& w) const { return (u128(u64((w * iv) >> 64) + 1) * mod) >> 64; }
 static constexpr u64 R2(u64 n) {
  u64 r= -u128(n) % n;
  return u128(r) * r % n;
 }
};
