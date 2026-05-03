#pragma once
#include "_common.hpp"
// 素朴な実装: 二分探索 with overflow safe pow。各クエリ O(64 · log A) 程度。
struct Solver {
 // overflow safe pow: a^b > limit なら true
 static bool pow_exceeds(u64 a, int b, u64 limit) {
  u64 r = 1;
  while (b) {
   if (b & 1) {
    if (a != 0 && r > limit / a) return true;
    r *= a;
    if (r > limit) return true;
   }
   b >>= 1;
   if (b) {
    if (a > 0 && a > limit / a) {
     // a*a > limit (loose bound). 続行するため a を limit+1 にしておく
     if (b) return true;
    }
    if (a > 0) a *= a;
   }
  }
  (void) r;
  return false;
 }
 // 真面目に: a^b の値を返す (overflow したら u64(-1))
 static u64 pow_safe(u64 a, int b) {
  u64 r = 1;
  for (int i = 0; i < b; ++i) {
   if (a != 0 && r > u64(-1) / a) return u64(-1);
   r *= a;
  }
  return r;
 }
 static u64 kth_root(u64 A, int k) {
  if (A == 0) return 0;
  if (k == 1) return A;
  // 二分探索: lo^k ≤ A < (lo+1)^k
  // ⌊(2^64-1)^{1/k}⌋ < 2^{ceil(64/k)} なので hi をその大きさに
  u64 hi = u64(1) << ((64 + k - 1) / k);
  u64 lo = 0;
  while (hi - lo > 1) {
   u64 mid = lo + (hi - lo) / 2;
   if (pow_safe(mid, k) <= A) lo = mid;
   else hi = mid;
  }
  return lo;
 }
 static std::string run(const std::string& input) {
  std::istringstream in(input);
  std::ostringstream out;
  int T; in >> T;
  while (T--) {
   u64 A;
   int k;
   in >> A >> k;
   out << kth_root(A, k) << '\n';
  }
  return std::move(out).str();
 }
};
