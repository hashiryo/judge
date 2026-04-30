#pragma once
#include "_common.hpp"
// 拡張 Euclid を繰り返して gcd を取る (Bezout 係数も計算するので普通の Euclid より重い)。
// 純粋な比較用、わざと遅い参照点。
//
// gcd 計算自体は u64 で行う (a >= 2^63 で i64 化が破綻するため)。
// Bezout 係数 x, y は i64 で追跡し、計算コストは載せる (オーバフローして
// 正しい値にならなくても gcd の正しさには影響しない、ベンチ目的)。
struct G {
 static inline u64 gcd(u64 a, u64 b) {
  i64 x= 1, y= 0, x1= 0, y1= 1;
  while (b) {
   u64 q= a / b;
   u64 t= a - q * b;
   a= b;
   b= t;
   i64 tx= x - (i64)q * x1;
   x= x1;
   x1= tx;
   i64 ty= y - (i64)q * y1;
   y= y1;
   y1= ty;
  }
  // x, y は使わないが計算コストを最適化で消されないよう sink する。
  asm volatile("" : : "r"(x), "r"(y));
  return a;
 }
};
