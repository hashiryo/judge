#pragma once
#include "_common.hpp"
// 試し割り (素朴版)。最大 1e18 だと 1e9 回試す必要があり TLE 必至だが、
// 比較対象として置いておく。
struct Factorize {
 static vector<u64> factor_one(u64 n) {
  vector<u64> fs;
  for (u64 p = 2; p * p <= n && p < 1000000; ++p) {
   while (n % p == 0) { fs.push_back(p); n /= p; }
  }
  // p^2 > n になるか p ≥ 1e6 まで来た場合: 残りは素数 1 つ (1e18 まで) または
  // 大きい合成数。素朴のままでは判定不能なので、ここでは「n > 1 なら 1 個の
  // 素数」として吐く (10^12 までは正しい、大きいケースは TLE/誤答 OK)。
  if (n > 1) fs.push_back(n);
  std::sort(fs.begin(), fs.end());
  return fs;
 }
 static vector<vector<u64>> run(const vector<u64>& qs) {
  vector<vector<u64>> ans;
  ans.reserve(qs.size());
  for (auto x : qs) ans.push_back(factor_one(x));
  return ans;
 }
};
