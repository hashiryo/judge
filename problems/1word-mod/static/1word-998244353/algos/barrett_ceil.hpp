#pragma once
#include "_common.hpp"
struct DIV {
 u32 b;
 constexpr DIV(u32 b): b(b), x(u32(-1) / b + 1) {}
 constexpr u32 mod(u32 a) const {
  u32 r= a - u64(u32(u64(a) * x >> 32)) * b;
  return r >> 31 ? r + b : r;
 }
private:
 u32 x;
};
