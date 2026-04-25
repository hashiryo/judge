#pragma once
// 2^20 < mod <= 2^41
struct MP {
    u64 mod;
    constexpr MP(u64 m): mod(m), x((u128(1) << 84) / m) {}
    constexpr u64 mul(u64 l, u64 r) const { return rem(u128(l) * r); }
    constexpr u64 set(u64 n) const { return n; }
    constexpr u64 get(u64 n) const { return n >= mod ? n - mod : n; }
    constexpr u64 plus(u64 l, u64 r) const { return l += r, l < (mod << 1) ? l : l - (mod << 1); }
    constexpr u64 diff(u64 l, u64 r) const { return l -= r, l >> 63 ? l + (mod << 1) : l; }
private:
    u64 x;
    constexpr u128 quo(u128 n) const { return (n * x) >> 84; }
    constexpr u64 rem(u128 n) const { return n - quo(n) * mod; }
};
