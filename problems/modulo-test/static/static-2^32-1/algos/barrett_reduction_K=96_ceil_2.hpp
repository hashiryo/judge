#pragma once
#include "_common.hpp"
struct MP {  // mod < 2^32
 u32 mod;
 constexpr MP(): mod(0), x(0) {}
 constexpr MP(u32 m): mod(m), x(((u128(1) << 96) - 1) / m + 1) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 static constexpr inline u32 set(u32 n) { return n; }
 static constexpr inline u32 get(u32 n) { return n; }
 static constexpr inline u32 norm(u32 n) { return n; }
 constexpr inline u32 plus(u64 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u64 l, u32 r) const { return l-= r, l >> 63 ? l + mod : l; }
private:
 u128 x;
 constexpr inline u32 quo(u64 n) const { return (x * n) >> 96; }
 constexpr inline u32 rem(u64 n) const { return n - u64(quo(n)) * mod; }
};
