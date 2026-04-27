#pragma once
#include "_common.hpp"
struct DIV {
 u32 b;
 constexpr DIV(u32 b): b(b), x(u64(-1) / b + 1) {}
 constexpr u32 mod(u32 a) const { return a - u64(u32(a * x >> 64)) * b; }
private:
 u64 x;
};
