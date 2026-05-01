#pragma once
#include "_common.hpp"
// std::bitset で行をパックして [A|I] のガウス・ジョルダン消去。可逆でなければ空 vector を返す。
struct Inv {
 static constexpr int MAXN = 4096;
 using BS = std::bitset<2 * MAXN>;
 static vector<string> run(int n, const vector<string>& a) {
  vector<BS> rows(n);
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) rows[i][j] = a[i][j] == '1';
   rows[i][n + i] = 1;
  }
  // forward elimination
  for (int j = 0; j < n; ++j) {
   int piv = -1;
   for (int i = j; i < n; ++i)
    if (rows[i][j]) { piv = i; break; }
   if (piv == -1) return {};
   if (piv != j) std::swap(rows[piv], rows[j]);
   for (int i = j + 1; i < n; ++i)
    if (rows[i][j]) rows[i] ^= rows[j];
  }
  // back substitution
  for (int j = n - 1; j >= 0; --j) {
   for (int i = 0; i < j; ++i)
    if (rows[i][j]) rows[i] ^= rows[j];
  }
  vector<string> out(n, string(n, '0'));
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j)
    if (rows[i][n + j]) out[i][j] = '1';
  return out;
 }
};
