#pragma once
#include "_common.hpp"
// mod < 2^32
struct MP {
    u32 mod;
    constexpr MP(u32 m): mod(m) {}
    constexpr u32 mul(u32 l, u32 r) const { return u64(l) * r % mod; }
    constexpr u32 set(u32 n) const { return n; }
    constexpr u32 get(u32 n) const { return n; }
    constexpr u32 plus(u64 l, u32 r) const { return l += r, l < mod ? l : l - mod; }
    constexpr u32 diff(u64 l, u32 r) const { return l -= r, l >> 63 ? l + mod : l; }
};
