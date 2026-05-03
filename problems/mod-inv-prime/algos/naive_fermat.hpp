#pragma once
#include "_common.hpp"
// Fermat per query: a^{-1} = a^{p-2} mod p
// p ≈ 10^9 で per-query ~30 mul (binary exp on 30-bit exponent)。最も素朴。
struct ModInv {
 static u32 mod_pow(u32 b, u32 e, u32 m) {
  u64 r = 1, x = b;
  while (e) { if (e & 1) r = r * x % m; x = x * x % m; e >>= 1; }
  return u32(r);
 }
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) {
   if (a == 0) ans.push_back(u32(-1));
   else ans.push_back(mod_pow(a, p - 2, p));
  }
  return ans;
 }
};
