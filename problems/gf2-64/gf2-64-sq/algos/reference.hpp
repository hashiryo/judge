#pragma once
#include "_common.hpp"
// 素朴な reference: bit-by-bit polynomial squaring + reduce。
//   sq(a) は GF(2) では各 bit i を bit 2i に spread (cross term ゼロ)。
//   reduction は bit-by-bit。
struct GF2_64Op {
 static u64 sq(u64 a) {
  u64 lo = 0, hi = 0;
  for (int i = 0; i < 64; ++i) {
   if ((a >> i) & 1) {
    int j = 2 * i;
    if (j < 64) lo ^= u64(1) << j;
    else hi ^= u64(1) << (j - 64);
   }
  }
  // reduce hi の各 bit を P(x) に従って xor
  for (int i = 63; i >= 0; --i) {
   if ((hi >> i) & 1) {
    hi ^= u64(1) << i;
    lo ^= IRRED_LOW << i;
    if (i > 0) hi ^= IRRED_LOW >> (64 - i);
   }
  }
  return lo;
 }
 static vector<u64> run(const vector<u64>& as) {
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = sq(as[i]);
  return ans;
 }
};
