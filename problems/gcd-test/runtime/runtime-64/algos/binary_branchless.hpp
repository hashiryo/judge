#pragma once
#include "_common.hpp"
// Stein を min/max で書いた branchless 風 GCD。
// gcc/clang は std::min/max を cmov に落とすので分岐除去が期待できる。
struct G {
 static inline u64 gcd(u64 a, u64 b) {
  if (!a) return b;
  if (!b) return a;
  int s= __builtin_ctzll(a | b);
  a>>= __builtin_ctzll(a);
  b>>= __builtin_ctzll(b);
  while (a != b) {
   u64 mn= std::min(a, b), mx= std::max(a, b);
   a= mn;
   b= (mx - mn) >> __builtin_ctzll(mx - mn);
  }
  return a << s;
 }
};
