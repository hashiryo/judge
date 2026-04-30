#pragma once
#include "_common.hpp"
// Library/mylib/number_theory/mod_inv.hpp 由来。
// Bezout 係数は 1 個 (x) しか追跡しない (mod 逆元計算には片方で十分のため)。
// 普通の extgcd の半分の追加計算量。
struct G {
 static inline u64 gcd(u64 a, u64 b) {
  i64 x= 1, y= 0, z= 0;
  u64 q= 0, c= 0;
  while (b) {
   q= a / b;
   z= x;
   x= y;
   y= z - y * (i64)q;
   c= a;
   a= b;
   b= c - b * q;
  }
  asm volatile("" : : "r"(x));
  return a;
 }
};
