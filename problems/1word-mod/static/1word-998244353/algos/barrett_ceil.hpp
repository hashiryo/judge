#pragma once
#include "_common.hpp"
struct DIV {
 static constexpr u32 X = u32(-1) / B + 1;
 constexpr DIV() {}
 constexpr u32 mod(u32 a) const {
  u32 r= a - u64(u32(u64(a) * X >> 32)) * B;
  return r >> 31 ? r + B : r;
 }
};
