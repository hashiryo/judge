#pragma once
#include "_common.hpp"
// Plantard 乗算 (32-bit、mod < 2^32/φ) で Fermat の小定理: a^{-1} = a^{p-2}
//
// p = 998244353 < 2^30 < 2^32/φ なので 32-bit Plantard が使える。
// 通常の % p と比べて 1 mul の reduce で済むため per-mul が高速。
//
// 参考: judge/problems/modulo-test/runtime/runtime-30/algos/plantard_2.hpp
//   と https://zenn.dev/yatyou/articles/f8a3969c6e9f7b
struct ModInv {
 struct MP {  // mod < 2^32/phi (φ = 黄金比)
  u32 mod;
  u32 r2;
  u64 iv;
  constexpr MP(u32 m): mod(m), r2(u32(-u128(m) % m)), iv(inv_(m)) {}
  static constexpr u64 inv_(u64 n, int e = 6, u64 x = 1) {
   return e ? inv_(n, e - 1, x * (2 - x * n)) : x;
  }
  constexpr u32 reduce(u64 w) const { return u32((u128((w * iv) | u32(-1)) * mod) >> 64); }
  constexpr u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
  constexpr u32 set(u32 n) const { return mul(n, r2); }
  constexpr u32 get(u32 n) const { return reduce(n); }
  constexpr u32 pow(u32 base, u32 e) const {
   u32 r = set(1);
   while (e) {
    if (e & 1) r = mul(r, base);
    base = mul(base, base);
    e >>= 1;
   }
   return r;
  }
 };
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  MP M(p);
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) {
   if (a == 0) { ans.push_back(u32(-1)); continue; }
   const u32 a_m = M.set(a);
   const u32 inv_m = M.pow(a_m, p - 2);
   ans.push_back(M.get(inv_m));
  }
  return ans;
 }
};
