#pragma once
#include "_common.hpp"
// std::bitset で拡張行列 [A|b] を RREF し、特解 + null space 基底を返す。
struct Solve {
 static constexpr int MAXM = 4096;
 using BS = std::bitset<MAXM + 1>;  // 拡張列込み
 static vector<string> run(int n, int m, const vector<string>& a, const vector<int>& bvec) {
  vector<BS> rows(n);
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < m; ++j) rows[i][j] = a[i][j] == '1';
   rows[i][m] = bvec[i] & 1;
  }
  // 前進消去 (RREF)
  vector<int> pivot_col(n, -1);  // pivot_col[i] = i 行目のピボット列 (なければ -1)
  int piv = 0;
  for (int j = 0; j < m && piv < n; ++j) {
   int sel = -1;
   for (int i = piv; i < n; ++i)
    if (rows[i][j]) { sel = i; break; }
   if (sel == -1) continue;
   if (sel != piv) std::swap(rows[piv], rows[sel]);
   for (int i = piv + 1; i < n; ++i)
    if (rows[i][j]) rows[i] ^= rows[piv];
   pivot_col[piv] = j;
   ++piv;
  }
  // 不整合チェック (0...0 | 1 行があれば解なし)
  for (int i = piv; i < n; ++i)
   if (rows[i][m]) return {};
  // 後退消去で reduced row echelon form に
  for (int i = piv - 1; i >= 0; --i) {
   int pc = pivot_col[i];
   for (int k = 0; k < i; ++k)
    if (rows[k][pc]) rows[k] ^= rows[i];
  }
  // pivot_col[i] = c なら x[c] = rows[i][m] - sum_{j: rows[i][j], j>c} x[j] (mod 2)
  // 自由変数 (pivot 列じゃない列) は free_vars に。
  vector<bool> is_pivot(m, false);
  for (int i = 0; i < piv; ++i) is_pivot[pivot_col[i]] = true;
  vector<int> free_vars;
  for (int j = 0; j < m; ++j)
   if (!is_pivot[j]) free_vars.push_back(j);
  int R = (int) free_vars.size();
  // 特解: x[pivot_col[i]] = rows[i][m]、自由変数は 0
  vector<string> out(R + 1, string(m, '0'));
  for (int i = 0; i < piv; ++i)
   if (rows[i][m]) out[0][pivot_col[i]] = '1';
  // 各 null space 基底: f 番目の自由変数 = 1、他 free=0、pivot 行から束縛変数を決定
  for (int b = 0; b < R; ++b) {
   int fc = free_vars[b];
   out[b + 1][fc] = '1';
   for (int i = 0; i < piv; ++i)
    if (rows[i][fc]) out[b + 1][pivot_col[i]] = '1';
  }
  return out;
 }
};
