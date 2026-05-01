#pragma once
#include "_common.hpp"
// 交代行列 A の Pfaffian を skew-symmetric Gaussian elimination で求める素朴版。
// 元提出と同じ手順 (列 i-1 の sub-diagonal pivot 探し → 行/列 swap →
//   奇数 i で res *= -A[i][i-1] → 行 j を行 i で消去) を行うが、AVX2 を
//   使わず u64 での add-mul で書く。O((2N)^3)。
struct Pfaffian {
 static u32 qpow_(u64 x, u64 y) {
  u64 r = 1;
  while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
  return (u32) r;
 }
 static u32 run(int N, const vector<vector<u32>>& M_in) {
  int n = 2 * N;
  vector<vector<u32>> A(n, vector<u32>(n));
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j) A[i][j] = M_in[i][j];
  u64 res = 1;
  for (int i = 1; i < n; ++i) {
   if (A[i][i - 1] == 0) {
    int sel = -1;
    for (int j = i + 1; j < n; ++j) if (A[j][i - 1]) { sel = j; break; }
    if (sel != -1) {
     std::swap(A[i], A[sel]);
     // 列 i, sel も swap (相似ではないが Pfaffian は行/列同時 swap で sign 反転)
     for (int k = i; k < n; ++k) std::swap(A[k][i], A[k][sel]);
     res = (MOD - res) % MOD;
    }
   }
   if (A[i][i - 1] == 0) return 0;
   if (i & 1) {
    // res *= -A[i][i-1]
    u32 v = MOD - A[i][i - 1];
    if (v == MOD) v = 0;
    res = res * v % MOD;
   }
   u32 inv = qpow_(A[i][i - 1], MOD - 2);
   for (int j = i + 1; j < n; ++j) {
    if (A[j][i - 1] == 0) continue;
    u32 f = (u32) ((u64) A[j][i - 1] * inv % MOD);
    u32 nf = (MOD - f) % MOD;
    for (int k = i - 1; k < n; ++k)
     A[j][k] = (u32) ((A[j][k] + (u64) nf * A[i][k]) % MOD);
   }
  }
  return (u32) res;
 }
};
