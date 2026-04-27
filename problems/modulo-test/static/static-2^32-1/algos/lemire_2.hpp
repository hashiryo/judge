#pragma once
#include "_common.hpp"
// https://lemire.me/blog/2019/02/08/faster-remainders-when-the-divisor-is-a-constant-beating-compilers-and-libdivide/
// https://min-25.hatenablog.com/entry/2017/08/20/171214
struct MP {  // mod < 2^32
 u32 mod;
 constexpr MP(): mod(0), s(0), x(0), mask(0) {}
 constexpr MP(u32 m): mod(m), s(31 - __builtin_clz(m - 1)), x((((u128(1) << 64) << s) - 1) / m + 1), mask(((u128(1) << 64) << s) - 1) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 static constexpr inline u32 set(u32 n) { return n; }
 static constexpr inline u32 get(u32 n) { return n; }
 static constexpr inline u32 norm(u32 n) { return n; }
 constexpr inline u32 plus(u64 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u64 l, u32 r) const { return l-= r, l >> 63 ? l + mod : l; }
private:
 u8 s;
 u64 x;
 u128 mask;
 constexpr inline u32 rem(u64 n) const { return (((u128(x) * n) & mask) * mod) >> 64 >> s; }
};
