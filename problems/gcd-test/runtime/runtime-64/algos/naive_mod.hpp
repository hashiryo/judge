#pragma once
#include "_common.hpp"
// 教科書的 Euclidean: a, b → a%b の繰り返し。
struct G {
 static inline u64 gcd(u64 a, u64 b) {
  while (b) {
   u64 t= a % b;
   a= b;
   b= t;
  }
  return a;
 }
};
