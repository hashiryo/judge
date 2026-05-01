#pragma once
#include "_common.hpp"
// 標準的な O(N^2 2^N) subset convolution:
//   - rank vectors A[i][k] = (popcount(i)==k) ? a[i] : 0 を作る (k=0..N)
//   - 各 k について SOS / OR-convolution の forward transform (Möbius zeta)
//   - 各 i について A_hat * B_hat を多項式積 (rank モード) で C_hat に詰める
//   - 各 k について逆変換 (Möbius)
//   - c[i] = C_hat[i][popcount(i)]
struct SubsetConv {
 static vector<u32> run(int N, const vector<u32>& a, const vector<u32>& b) {
  int sz = 1 << N;
  vector<vector<u64>> A(N + 1, vector<u64>(sz, 0));
  vector<vector<u64>> B(N + 1, vector<u64>(sz, 0));
  for (int i = 0; i < sz; ++i) {
   int p = __builtin_popcount(i);
   A[p][i] = a[i];
   B[p][i] = b[i];
  }
  auto zeta = [&](vector<u64>& f) {
   for (int i = 0; i < N; ++i)
    for (int j = 0; j < sz; ++j)
     if (j & (1 << i)) f[j] = (f[j] + f[j ^ (1 << i)]) % MOD;
  };
  auto mobius = [&](vector<u64>& f) {
   for (int i = 0; i < N; ++i)
    for (int j = 0; j < sz; ++j)
     if (j & (1 << i)) f[j] = (f[j] + MOD - f[j ^ (1 << i)]) % MOD;
  };
  for (int k = 0; k <= N; ++k) zeta(A[k]);
  for (int k = 0; k <= N; ++k) zeta(B[k]);
  vector<vector<u64>> C(N + 1, vector<u64>(sz, 0));
  for (int i = 0; i < sz; ++i) {
   for (int p = 0; p <= N; ++p) {
    u64 sum = 0;
    for (int q = 0; q <= p; ++q)
     sum = (sum + A[q][i] * B[p - q][i]) % MOD;
    C[p][i] = sum;
   }
  }
  for (int k = 0; k <= N; ++k) mobius(C[k]);
  vector<u32> c(sz);
  for (int i = 0; i < sz; ++i)
   c[i] = (u32) C[__builtin_popcount(i)][i];
  return c;
 }
};
