#pragma once
#include "_common.hpp"
struct MP {  // mod < 2^30
 u32 mod;
 constexpr MP(): mod(0), s(0), mo(0), r2(0), mask(0), iv1(0), iv2(0) {}
 MP(u32 m): mod(m), s(__builtin_ctz(m)), mo(m >> s), r2(-u64(mo) % mo), mask((1ull << (s + 32)) - 1) {
  u64 iv= inv(mo);
  iv1= iv * -u64(u32(-1)), iv2= iv * -u64(u32(-r2));
 }
 constexpr inline u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return transform(n); }
 constexpr inline u32 get(u32 n) const { return n= reduce(n), n >= mod ? n - mod : n; }
 constexpr inline u32 norm(u32 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < (mod << 1) ? l : l - (mod << 1); }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> (32 - 1) ? l + (mod << 1) : l; }
private:
 u8 s;
 u32 mo, r2;
 u64 mask, iv1, iv2;
 static constexpr u64 inv(u32 n, int e= 6, u64 x= 1) { return e ? inv(n, e - 1, x * (2 - x * n)) : x; }
 constexpr inline u32 reduce(u64 w) const { return u32(w >> 32) + mod - u32((((w * iv1) & mask) * mo) >> 32); }
 constexpr inline u32 transform(u32 n) const { return u32(u64(n) * r2 >> 32) + mod - u32((((n * iv2) & mask) * mo) >> 32); }
};
