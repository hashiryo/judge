#pragma once
#include "_common.hpp"
struct DIV {
 u64 b;
 constexpr DIV(u64 b): b(b), x(u128(-1) / b + 1) {}
 constexpr u64 mod(u64 a) const { return a - u128(u64(a * x >> 64)) * b; }
private:
 u128 x;
};
