#pragma once
#include "_common.hpp"
// https://zenn.dev/yatyou/articles/f8a3969c6e9f7b
struct MP {  // mod < 2^32/phi
 u32 mod;
 constexpr MP(): mod(0), iv(0), r2(0) {}
 constexpr MP(u32 m): mod(m), iv(inv(m)), r2(-u128(mod) % mod) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return mul(n, r2); }
 constexpr inline u32 get(u32 n) const { return n= reduce(n), n >= mod ? n - mod : n; }
 constexpr inline u32 norm(u32 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + mod : l; }
private:
 u64 iv, r2;
 static constexpr u64 inv(u64 n, int e= 6, u64 x= 1) { return e ? inv(n, e - 1, x * (2 - x * n)) : x; }
 constexpr inline u32 reduce(u64 w) const { return (u64(u32((w * iv) >> 32) + 1) * mod) >> 32; }
};
