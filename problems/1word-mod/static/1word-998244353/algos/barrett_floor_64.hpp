#pragma once
#include "_common.hpp"
struct DIV {
 static constexpr u64 X = -u64(B) / B + 1;
 constexpr DIV() {}
 constexpr u32 mod(u32 a) const {
  u32 r= a - u32(u128(a) * X >> 64) * B;
  return r >= B ? r - B : r;
 }
};
