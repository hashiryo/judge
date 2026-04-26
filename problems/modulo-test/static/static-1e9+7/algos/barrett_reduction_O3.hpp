#pragma once
#include "_common.hpp"
#pragma GCC optimize("O3")
struct MP {  // 2^20 < mod <= 2^41
 u64 mod;
 constexpr MP(): mod(0), x(0) {}
 constexpr MP(u64 m): mod(m), x((u128(1) << 84) / m) {}
 constexpr inline u64 mul(u64 l, u64 r) const { return rem(u128(l) * r); }
 static constexpr inline u64 set(u64 n) { return n; }
 constexpr inline u64 get(u64 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u64 norm(u64 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u64 plus(u64 l, u64 r) const { return l+= r, l < (mod << 1) ? l : l - (mod << 1); }
 constexpr inline u64 diff(u64 l, u64 r) const { return l-= r, l >> 63 ? l + (mod << 1) : l; }
private:
 u64 x;
 constexpr inline u128 quo(const u128& n) const { return (n * x) >> 84; }
 constexpr inline u64 rem(const u128& n) const { return n - quo(n) * mod; }
};
