#pragma once
#include "_common.hpp"
// 拡張行列 [A | I] を Gauss-Jordan で [I | A^{-1}] にする素朴版。
struct Inverse {
 static InverseResult run(int n, const vector<vector<u32>>& a_in) {
  auto qpow = [](u64 x, u64 y) {
   u64 r = 1;
   while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
   return (u32) r;
  };
  vector<vector<u32>> M_(n, vector<u32>(2 * n, 0));
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) M_[i][j] = a_in[i][j];
   M_[i][n + i] = 1;
  }
  for (int i = 0; i < n; ++i) {
   int sel = -1;
   for (int j = i; j < n; ++j) if (M_[j][i]) { sel = j; break; }
   if (sel == -1) return {false, {}};
   if (sel != i) std::swap(M_[i], M_[sel]);
   u32 inv = qpow(M_[i][i], MOD - 2);
   for (int k = 0; k < 2 * n; ++k)
    M_[i][k] = (u32) ((u64) M_[i][k] * inv % MOD);
   for (int j = 0; j < n; ++j) {
    if (j == i || !M_[j][i]) continue;
    u32 f = MOD - M_[j][i];
    for (int k = 0; k < 2 * n; ++k)
     M_[j][k] = (u32) ((M_[j][k] + (u64) f * M_[i][k]) % MOD);
   }
  }
  vector<vector<u32>> ans(n, vector<u32>(n));
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j) ans[i][j] = M_[i][n + j];
  return {true, std::move(ans)};
 }
};
