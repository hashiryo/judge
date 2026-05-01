#pragma once
#include "_common.hpp"
// std::bitset で行をパックして GF(2) ガウス消去。pivot が見つからなければ det=0。
struct Det {
 static constexpr int MAXN = 4096;
 using BS = std::bitset<MAXN>;
 static int run(int n, const vector<string>& a) {
  vector<BS> rows(n);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j) rows[i][j] = a[i][j] == '1';
  for (int j = 0; j < n; ++j) {
   int piv = -1;
   for (int i = j; i < n; ++i)
    if (rows[i][j]) { piv = i; break; }
   if (piv == -1) return 0;
   if (piv != j) std::swap(rows[piv], rows[j]);
   for (int i = j + 1; i < n; ++i)
    if (rows[i][j]) rows[i] ^= rows[j];
  }
  return 1;
 }
};
