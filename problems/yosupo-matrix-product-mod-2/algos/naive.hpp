#pragma once
#include "_common.hpp"
// std::bitset で行をパックして行列積。
// c[i][j] = parity of (row_i_of_A AND col_j_of_B)。
// B を「列ごとに bitset」で持って AND→popcount で計算する。
struct Mul {
 static constexpr int MAXN = 4096;
 using BS = std::bitset<MAXN>;
 static vector<string> run(int n, int m, int k, const vector<string>& a, const vector<string>& b) {
  vector<BS> ar(n), bc(k);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < m; ++j) ar[i][j] = a[i][j] == '1';
  for (int i = 0; i < m; ++i)
   for (int j = 0; j < k; ++j) bc[j][i] = b[i][j] == '1';
  vector<string> res(n, string(k, '0'));
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < k; ++j) {
    res[i][j] = ((ar[i] & bc[j]).count() & 1) ? '1' : '0';
   }
  }
  return res;
 }
};
