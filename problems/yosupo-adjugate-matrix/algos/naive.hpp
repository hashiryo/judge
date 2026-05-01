#pragma once
#include "_common.hpp"
// adj(A) = det(A) * A^{-1} (det != 0 のとき)。素朴に Gauss-Jordan で逆行列と
// det を同時に求める。det == 0 の場合は cp-algo の hack と同じく
// 「(n+1)×(n+1) に拡張して逆行列の (n,n) cofactor から復元」する方式で計算
// する手もあるが、素朴版では det != 0 の場合のみ A^{-1} 経路、それ以外は
// 余因子展開を直接やる O(N^4) フォールバック。
struct Adjugate {
 static u32 qpow_(u64 x, u64 y) {
  u64 r = 1;
  while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
  return (u32) r;
 }
 // det を求めつつ消去 (素朴 Gauss)
 static u32 det_of(vector<vector<u32>> a, int n) {
  bool neg = false;
  for (int i = 0; i < n; ++i) {
   int sel = -1;
   for (int j = i; j < n; ++j) if (a[j][i]) { sel = j; break; }
   if (sel == -1) return 0;
   if (sel != i) { std::swap(a[i], a[sel]); neg = !neg; }
   u32 inv = qpow_(a[i][i], MOD - 2);
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
 static vector<vector<u32>> run(int n, const vector<vector<u32>>& a_in) {
  vector<vector<u32>> ans(n, vector<u32>(n, 0));
  // det(A) を計算
  u32 D = det_of(a_in, n);
  if (D != 0) {
   // [A | I] → [I | A^{-1}], adj = det * A^{-1}
   vector<vector<u32>> M_(n, vector<u32>(2 * n, 0));
   for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) M_[i][j] = a_in[i][j];
    M_[i][n + i] = 1;
   }
   for (int i = 0; i < n; ++i) {
    int sel = -1;
    for (int j = i; j < n; ++j) if (M_[j][i]) { sel = j; break; }
    if (sel != i) std::swap(M_[i], M_[sel]);
    u32 inv = qpow_(M_[i][i], MOD - 2);
    for (int k = 0; k < 2 * n; ++k) M_[i][k] = (u32) ((u64) M_[i][k] * inv % MOD);
    for (int j = 0; j < n; ++j) {
     if (j == i || !M_[j][i]) continue;
     u32 f = MOD - M_[j][i];
     for (int k = 0; k < 2 * n; ++k)
      M_[j][k] = (u32) ((M_[j][k] + (u64) f * M_[i][k]) % MOD);
    }
   }
   for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) ans[i][j] = (u32) ((u64) D * M_[i][n + j] % MOD);
   return ans;
  }
  // det == 0: 各 (i,j) の余因子 (-1)^{i+j} M_{ji} を直接計算 (O(N^4))。
  // adj(A)[i][j] = (-1)^{i+j} det(A の i 列 j 行を除いた (n-1)×(n-1) 行列)
  // ※入力では (i,j) 表記が転置されるので注意: ans[i][j] = (-1)^{i+j} M_{ji}
  if (n == 1) { ans[0][0] = 1; return ans; }
  vector<vector<u32>> minor_(n - 1, vector<u32>(n - 1));
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) {
    // M_{ji} = det(A から j 行目, i 列目を除いた行列)
    int rr = 0;
    for (int r = 0; r < n; ++r) {
     if (r == j) continue;
     int cc = 0;
     for (int c = 0; c < n; ++c) {
      if (c == i) continue;
      minor_[rr][cc] = a_in[r][c];
      ++cc;
     }
     ++rr;
    }
    u32 m = det_of(minor_, n - 1);
    if ((i + j) & 1) m = m ? MOD - m : 0;
    ans[i][j] = m;
   }
  }
  return ans;
 }
};
