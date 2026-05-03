#pragma once
#include "_common.hpp"
// 拡張ユークリッド per query: a · x + p · y = 1 を解いて x = a^{-1}
// 期待 ~ log φ p ≈ 30 iter、各 iter で 1 div + 数 sub。Fermat 比 ~3-5× 速い見込み。
struct ModInv {
 static u32 inv_extgcd(u32 a, u32 m) {
  if (a == 0) return u32(-1);
  u32 ori_m = m;
  i32 b = 1, c = 0;
  while (a != 1) {
   if (a == 0) return u32(-1);
   u32 d = m / a;
   c -= b * (i32) d;
   m -= a * d;
   std::swap(a, m); std::swap(b, c);
  }
  return b < 0 ? u32(b + (i32) ori_m) : u32(b);
 }
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) ans.push_back(inv_extgcd(a, p));
  return ans;
 }
};
