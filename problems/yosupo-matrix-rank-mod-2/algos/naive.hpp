#pragma once
#include "_common.hpp"
// std::bitset で行をパックして GF(2) ガウス消去。pivot 行数 = ランク。
struct Rank {
 static constexpr int MAXM = 4096;
 using BS = std::bitset<MAXM>;
 static int run(int n, int m, const vector<string>& a) {
  vector<BS> rows(n);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < m; ++j) rows[i][j] = a[i][j] == '1';
  int piv = 0;
  for (int j = 0; j < m && piv < n; ++j) {
   int sel = -1;
   for (int i = piv; i < n; ++i)
    if (rows[i][j]) { sel = i; break; }
   if (sel == -1) continue;
   if (sel != piv) std::swap(rows[piv], rows[sel]);
   for (int i = piv + 1; i < n; ++i)
    if (rows[i][j]) rows[i] ^= rows[piv];
   ++piv;
  }
  return piv;
 }
};
