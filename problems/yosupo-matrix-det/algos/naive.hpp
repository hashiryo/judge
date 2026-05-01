#pragma once
#include "_common.hpp"
// mod 998244353 (素数) 上の素朴 Gauss 消去で det を計算。
struct Det {
 static u32 run(int n, const vector<vector<u32>>& a_in) {
  vector<vector<u32>> a = a_in;
  auto qpow = [](u64 x, u64 y) {
   u64 r = 1;
   while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
   return (u32) r;
  };
  bool neg = false;
  for (int i = 0; i < n; ++i) {
   int sel = -1;
   for (int j = i; j < n; ++j)
    if (a[j][i]) { sel = j; break; }
   if (sel == -1) return 0;
   if (sel != i) { std::swap(a[i], a[sel]); neg = !neg; }
   u32 inv = qpow(a[i][i], MOD - 2);
   for (int j = i + 1; j < n; ++j) {
    u32 f = (u32) ((u64) a[j][i] * inv % MOD);
    if (!f) continue;
    f = MOD - f;
    for (int k = i; k < n; ++k)
     a[j][k] = (u32) ((a[j][k] + (u64) f * a[i][k]) % MOD);
   }
  }
  u64 res = neg ? MOD - 1 : 1;
  for (int i = 0; i < n; ++i) res = res * a[i][i] % MOD;
  return (u32) res;
 }
};
