#pragma once
#include "_common.hpp"
// Library/mylib/number_theory/binary_gcd.hpp の写し。
// branchless: swap せず三項演算子 (cmov) で進める。両方を最初に odd 化、ループ内で
// d = a - b の bsf で d をシフトしてから次のイテへ。
struct G {
 static inline u64 gcd(u64 a, u64 b) {
  if (a == 0 || b == 0) return a + b;
  int n= __builtin_ctzll(a), m= __builtin_ctzll(b), s= 0;
  for (a>>= n, b>>= m; a != b;) {
   u64 d= a - b;
   bool f= a > b;
   s= __builtin_ctzll(d);
   b= f ? b : a;
   a= (f ? d : -d) >> s;
  }
  return a << std::min(n, m);
 }
};
