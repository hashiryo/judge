#pragma once
#include "_common.hpp"
// 拡張ユークリッド互除法の符号無し整数版。
//
// 通常の ext_gcd は (a, b) を符号付きで持つが、以下の漸化式で常に非負に保てる:
//   不変条件: a · n ≡ x  (mod m),  b · n ≡ -y  (mod m)  (a, b ∈ [0, m))
//   1 行目に 2 行目の t 倍を加える操作 (B): a += b · t, y %= x
//   2 行目に 1 行目の t 倍を加える操作 (A): b += a · t, x %= y
//   x = 1 なら n^{-1} ≡ a, y = 1 なら n^{-1} ≡ -b ≡ m - b
//
// 参考: https://zenn.dev/mizar/articles/ea1aaa320a4b9d
//   (「一般の m に関するモジュラ逆数 (変則・拡張ユークリッド互除法)」)
struct ModInv {
 static u32 inv_unsigned(u32 n, u32 m) {
  if (m == 1) return 0;
  if (n == 0) return u32(-1);
  u64 a = 1, b = 0;
  u32 x = n, y = m;
  while (true) {
   if (x == 0) return u32(-1);  // gcd != 1
   if (x == 1) return u32(a % m);
   b += a * (y / x);
   y %= x;
   if (y == 0) return u32(-1);
   if (y == 1) return u32(m - b % m);
   a += b * (x / y);
   x %= y;
  }
 }
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) ans.push_back(inv_unsigned(a, p));
  return ans;
 }
};
