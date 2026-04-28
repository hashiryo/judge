#pragma once
#include "_common.hpp"
template <class Uint> constexpr inline Uint mod_inv(Uint a, Uint mod) {
 std::make_signed_t<Uint> x= 1, y= 0, z= 0;
 for(Uint q= 0, b= mod, c= 0; b;) z= x, x= y, y= z - y * (q= a / b), c= a, a= b, b= c - b * q;
 return assert(a == 1), x < 0 ? mod - (-x) % mod : x % mod;
}
struct MP {  // mod < 2^30
 u32 mod;
 constexpr MP(): mod(0), s(0), mo(0), d(0), mask(0), c(0), iv(0) {}
 MP(u32 m): mod(m), s(__builtin_ctz(m)), mo(m >> s), d(mod_inv(1u << s, mo)), mask((1u << s) - 1), c(((((1ull << (64 - s)) - d) % mo) << s) + 1), iv(inv(mo)) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return mul(n, c); }
 constexpr inline u32 get(u32 n) const { return n= reduce(n), n >= mod ? n - mod : n; }
 constexpr inline u32 norm(u32 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < (mod << 1) ? l : l - (mod << 1); }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + (mod << 1) : l; }
private:
 u8 s;
 u32 mo, d, mask, c, iv;
 static constexpr u32 inv(u32 n, int e= 5, u32 x= 1) { return e ? inv(n, e - 1, x * (2 - x * n)) : x; }
 constexpr inline u32 reduce(u64 w) const {
  u32 p= u32(w) & mask;
  u64 x= (w >> s) + u64(p) * d;
  u32 y= p << (32 - s);
  return (x + u64(mo) * (iv * (y - u32(x)))) >> (32 - s);
 }
};
