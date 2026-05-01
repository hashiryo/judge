#pragma once
#include "_common.hpp"
// mod 998244353 上の Ax=b を素朴 Gauss 消去で解く。
// 拡張行列 [A|b] を RREF まで持っていき、
//   - 主元が b 列 (右端) にある行があれば解なし
//   - そうでなければ自由変数を 0 とした特解 sol、自由変数 ↔ basis 行を構築
struct Solver {
 static SolveResult run(int n, int m, const vector<vector<u32>>& A_in, const vector<u32>& b_in) {
  auto qpow = [](u64 x, u64 y) {
   u64 r = 1;
   while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
   return (u32) r;
  };
  // [A | b] の拡張行列を作る (n × (m+1))
  vector<vector<u32>> M_(n, vector<u32>(m + 1, 0));
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < m; ++j) M_[i][j] = A_in[i][j];
   M_[i][m] = b_in[i];
  }
  vector<int> pivot_col(n, -1);
  vector<int> col_pivot(m + 1, -1);
  int row = 0;
  for (int col = 0; col < m && row < n; ++col) {
   int sel = -1;
   for (int j = row; j < n; ++j)
    if (M_[j][col]) { sel = j; break; }
   if (sel == -1) continue;
   if (sel != row) std::swap(M_[row], M_[sel]);
   u32 inv = qpow(M_[row][col], MOD - 2);
   // 行を 1/pivot 倍
   for (int k = col; k <= m; ++k)
    M_[row][k] = (u32) ((u64) M_[row][k] * inv % MOD);
   // 他行から消去 (RREF)
   for (int j = 0; j < n; ++j) {
    if (j == row || !M_[j][col]) continue;
    u32 f = MOD - M_[j][col];
    for (int k = col; k <= m; ++k)
     M_[j][k] = (u32) ((M_[j][k] + (u64) f * M_[row][k]) % MOD);
   }
   pivot_col[row] = col;
   col_pivot[col] = row;
   ++row;
  }
  // 解なし判定: 主元が m 列 (b 側) のみの行
  for (int j = 0; j < n; ++j) {
   bool zero_lhs = true;
   for (int k = 0; k < m; ++k) if (M_[j][k]) { zero_lhs = false; break; }
   if (zero_lhs && M_[j][m]) return {false, {}, {}};
  }
  // 特解: 自由変数 = 0, ピボット列 i = M_[row(i)][m]
  vector<u32> sol(m, 0);
  for (int r = 0; r < row; ++r) {
   int c = pivot_col[r];
   if (c >= 0) sol[c] = M_[r][m];
  }
  // basis: 自由変数ごとに 1 つ basis vector
  vector<vector<u32>> basis;
  for (int c = 0; c < m; ++c) {
   if (col_pivot[c] != -1) continue;
   vector<u32> v(m, 0);
   v[c] = 1;
   for (int r = 0; r < row; ++r) {
    int pc = pivot_col[r];
    if (pc < 0) continue;
    // x_{pc} + ... + M_[r][c] * x_c + ... = 0  → x_{pc} = -M_[r][c]
    if (M_[r][c]) v[pc] = MOD - M_[r][c];
   }
   basis.push_back(std::move(v));
  }
  return {true, std::move(sol), std::move(basis)};
 }
};
