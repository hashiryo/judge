#pragma once
#include "_common.hpp"
// Stein のアルゴリズム (binary GCD)。除算を使わず ctz と引き算で進める。
// gcd(2a, 2b) = 2 gcd(a,b), gcd(2a, b) = gcd(a, b) (b 奇), gcd(a, b) = gcd(|a-b|, min(a,b)) (両奇)。
struct G {
 static inline u64 gcd(u64 a, u64 b) {
  if (!a) return b;
  if (!b) return a;
  int s= __builtin_ctzll(a | b);
  a>>= __builtin_ctzll(a);
  do {
   b>>= __builtin_ctzll(b);
   if (a > b) std::swap(a, b);
   b-= a;
  } while (b);
  return a << s;
 }
};
