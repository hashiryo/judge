#pragma once
#include "_common.hpp"
// 素朴な実装: 二分探索 with overflow safe pow。各クエリ O(64 · log A) 程度。
struct Solver {
 // a^b ≤ A かどうかを overflow safe に判定
 static bool pow_le(u64 a, int b, u64 A) {
  u64 r = 1;
  for (int i = 0; i < b; ++i) {
   if (a != 0 && r > A / a) return false;  // overflow しそう or a^b > A
   r *= a;
   if (r > A) return false;
  }
  return true;
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
   if (pow_le(mid, k, A)) lo = mid;
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
