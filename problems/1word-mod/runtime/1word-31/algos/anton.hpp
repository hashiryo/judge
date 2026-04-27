#pragma once
#include "_common.hpp"
struct DIV {
 u32 b;
 constexpr DIV(u32 b): b(b), x(u64(-1) / b + 1) {}
 constexpr u32 mod(u32 a) const { return u64(u32(x * a >> 32) + 1) * b >> 32; }
private:
 u64 x;
};
