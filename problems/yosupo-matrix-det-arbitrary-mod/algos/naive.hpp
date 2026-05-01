#pragma once
#include "_common.hpp"
// 任意 mod 用 GF gauss: 体ではないので「ピボットを互除法で 0 にする」方式。
// 行 j の i 列目を、行 i の i 列目を法として消す: 互除法 (繰り返し除算 + swap)。
// 教科書通り、O(N^3 log M)。
struct Det {
 static u32 run(int n, u32 mod, const vector<vector<u32>>& a_in) {
  if (mod == 1) return 0;
  vector<vector<u32>> a = a_in;
  u32 res = 1;
  bool neg = false;
  for (int i = 0; i < n; ++i) {
   for (int j = i + 1; j < n; ++j) {
    while (a[i][i] != 0) {
     u32 q = a[j][i] / a[i][i];
     u32 inv = (u32) (mod - q);
     if (inv != mod) {
      for (int x = i; x < n; ++x) a[j][x] = (u32) ((a[j][x] + (u64) inv * a[i][x]) % mod);
     }
     std::swap(a[i], a[j]);
     neg = !neg;
    }
    std::swap(a[i], a[j]);
    neg = !neg;
   }
  }
  for (int i = 0; i < n; ++i) res = (u32) ((u64) res * a[i][i] % mod);
  if (neg && res != 0) res = mod - res;
  return res;
 }
};
