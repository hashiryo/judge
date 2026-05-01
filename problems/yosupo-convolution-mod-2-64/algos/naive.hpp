#pragma once
#include "_common.hpp"
// O(N*M) 素朴畳み込み。u64 のラップアラウンドで自動的に mod 2^64。
struct Conv {
 static vector<u64> run(const vector<u64>& a, const vector<u64>& b) {
  int n= (int)a.size(), m= (int)b.size();
  if (!n || !m) return {};
  vector<u64> c(n + m - 1, 0);
  for (int i= 0; i < n; ++i)
   for (int j= 0; j < m; ++j) c[i + j]+= a[i] * b[j];
  return c;
 }
};
