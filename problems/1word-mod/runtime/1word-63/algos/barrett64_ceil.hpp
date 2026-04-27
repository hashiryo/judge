#pragma once
#include "_common.hpp"
struct DIV {
 u64 b;
 constexpr DIV(u64 b): b(b), x(u64(-1) / b + 1) {}
 constexpr u64 mod(u64 a) const {
  u64 r= a - u128(u64(u128(a) * x >> 32)) * b;
  return r >> 31 ? r + b : r;
 }
private:
 u64 x;
};
