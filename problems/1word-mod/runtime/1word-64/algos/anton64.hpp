#pragma once
#include "_common.hpp"
struct DIV {
 u64 b;
 constexpr DIV(u64 b): b(b), x(u128(-1) / b + 1) {}
 constexpr u64 mod(u64 a) const { return u128(u64(x * a >> 32) + 1) * b >> 32; }
private:
 u128 x;
};
