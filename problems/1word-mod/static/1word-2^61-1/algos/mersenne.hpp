#pragma once
#include "_common.hpp"
// B = 2^61 - 1 の Mersenne 専用: a < 2^63 のとき t = (a & B) + (a >> 61) は t < 2 * B。
struct DIV {
 constexpr DIV() {}
 constexpr u64 mod(u64 a) const {
  u64 t= (a & B) + (a >> 61);
  return t >= B ? t - B : t;
 }
};
