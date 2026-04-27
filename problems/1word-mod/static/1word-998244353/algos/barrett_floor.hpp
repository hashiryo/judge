#pragma once
#include "_common.hpp"
struct DIV {
 static constexpr u32 X = -u32(B) / B + 1;
 constexpr DIV() {}
 constexpr u32 mod(u32 a) const {
  u32 r= a - u64(u32(u64(a) * X >> 32)) * B;
  return r >= B ? r - B : r;
 }
};
