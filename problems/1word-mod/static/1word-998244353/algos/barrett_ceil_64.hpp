#pragma once
#include "_common.hpp"
struct DIV {
 static constexpr u64 X = u64(-1) / B + 1;
 constexpr DIV() {}
 constexpr u32 mod(u32 a) const { return a - u32(u128(a) * X >> 64) * B; }
};
