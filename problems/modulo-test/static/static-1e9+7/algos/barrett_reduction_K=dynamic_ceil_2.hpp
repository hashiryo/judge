#pragma once
#include "_common.hpp"
struct MP {  // mod < 2^31
 u32 mod;
 constexpr MP(): mod(0), s(0), x(0) {}
 constexpr MP(u32 m): mod(m), s(31 - __builtin_clz(m - 1)), x((((u128(1) << 64) << s) - 1) / m + 1) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 static constexpr inline u32 set(u32 n) { return n; }
 static constexpr inline u32 get(u32 n) { return n; }
 static constexpr inline u32 norm(u32 n) { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + mod : l; }
private:
 u8 s;
 u64 x;
 constexpr inline u32 quo(u64 n) const { return (u128(x) * n) >> 64 >> s; }
 constexpr inline u32 rem(u64 n) const { return n - u64(quo(n)) * mod; }
};
