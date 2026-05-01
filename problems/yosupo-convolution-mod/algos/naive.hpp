#pragma once
#include "_common.hpp"
// O(N*M) 素朴畳み込み。サンプル/小ケース動作確認用。大きい N では TLE するので
// 比較ラインとしての参照点。
struct Conv {
 static vector<u32> run(const vector<u32>& a, const vector<u32>& b) {
  int n= (int)a.size(), m= (int)b.size();
  if (!n || !m) return {};
  vector<u32> c(n + m - 1, 0);
  for (int i= 0; i < n; ++i)
   for (int j= 0; j < m; ++j) c[i + j]= (u32)((c[i + j] + (u64)a[i] * b[j]) % MOD);
  return c;
 }
};
