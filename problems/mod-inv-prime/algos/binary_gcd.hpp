#pragma once
#include "_common.hpp"
// Binary GCD ベースのモジュラ逆数。
//
// 不変条件: a · n ≡ x (mod m),  b · n ≡ y (mod m)  (a, b ∈ [0, m))
// 操作:
//   x が偶数: x /= 2、a を mod m で半分にする (a 偶 → a/2、a 奇 → (a+m)/2)
//   y が偶数: 同様に y /= 2、b を mod m で半分
//   x ≤ y:    y -= x, b ← b - a (mod m)
//   x > y:    swap して上を実施
//
// y が 0 になったとき、x = gcd(n, m)。
// 素数 m なら gcd = 1 なので a = n^{-1}。
//
// 大きい商に対しては割り算 1 回で多倍進めるユークリッドが速いが、本実装では
// 常に 1 単位ずつ減らすので比例して遅い。ただしハードウェア div を使わない
// pure shift+sub の構造なので低レイテンシ命令だけで構成される。
//
// 参考: https://zenn.dev/mizar/articles/ea1aaa320a4b9d
//   (「式変形手順の具体例 ( m が 33 以上の奇数 の場合、バイナリGCD法の複合 )」
//    の単純化版 — A.B. ルールを使わず C.D. のみ)
struct ModInv {
 static u32 inv_bgcd(u32 n, u32 m) {
  if (n == 0) return u32(-1);
  if (n == 1) return 1;
  u32 x = n, y = m;
  u32 a = 1, b = 0;
  // x の偶数因子を抜く (a も mod m で同じ回数 halving)
  while ((x & 1) == 0) {
   x >>= 1;
   a = (a & 1) ? (u32) ((u64(a) + m) >> 1) : (a >> 1);
  }
  while (y != 0) {
   // y を奇数化
   while ((y & 1) == 0) {
    y >>= 1;
    b = (b & 1) ? (u32) ((u64(b) + m) >> 1) : (b >> 1);
   }
   // x, y 共に奇数にして大小比較・減算
   if (x > y) { std::swap(x, y); std::swap(a, b); }
   y -= x;
   b = (b >= a) ? (b - a) : (b + m - a);
  }
  return x == 1 ? a : u32(-1);
 }
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) ans.push_back(inv_bgcd(a, p));
  return ans;
 }
};
