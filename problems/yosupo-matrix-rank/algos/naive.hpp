#pragma once
#include "_common.hpp"
// mod 998244353 (素数) 上の素朴 Gauss 消去で rank を計算。
// N < M なら転置すると pivot 列探索が短縮されるが、行サイズも変わるので
// ここでは特別扱いせず素朴に行く。
struct Rank {
 static int run(int n, int m, const vector<vector<u32>>& a_in) {
  vector<vector<u32>> a = a_in;
  auto qpow = [](u64 x, u64 y) {
   u64 r = 1;
   while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
   return (u32) r;
  };
  int rank = 0;
  int col = 0;
  for (int row = 0; row < n && col < m; ) {
   int sel = -1;
   for (int j = row; j < n; ++j)
    if (a[j][col]) { sel = j; break; }
   if (sel == -1) { ++col; continue; }
   if (sel != row) std::swap(a[row], a[sel]);
   u32 inv = qpow(a[row][col], MOD - 2);
   for (int j = 0; j < n; ++j) {
    if (j == row || !a[j][col]) continue;
    u32 f = MOD - (u32) ((u64) a[j][col] * inv % MOD);
    for (int k = col; k < m; ++k)
     a[j][k] = (u32) ((a[j][k] + (u64) f * a[row][k]) % MOD);
   }
   ++rank; ++row; ++col;
  }
  return rank;
 }
};
