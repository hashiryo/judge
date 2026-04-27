#pragma once
#include "_common.hpp"
struct DIV {
 u32 b;
 constexpr DIV(u32 b): b(b), x(u64(-1) / b + 1) {}
 constexpr u32 mod(u32 a) const { return u128(u64(x * a)) * b >> 64; }
private:
 u64 x;
};
