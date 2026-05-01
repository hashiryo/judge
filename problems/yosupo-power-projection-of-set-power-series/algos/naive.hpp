#pragma once
#include "_common.hpp"
// 素朴: b = a^m を毎ステップ subset_convolution で更新しながら w·b を集める。
// O(M * N^2 * 2^N)。M ≤ 1e5, N ≤ 20 なら大きいケースでは TLE する。
struct PowerProj {
 static vector<u32> conv_(int N, const vector<u32>& a, const vector<u32>& b) {
  int sz = 1 << N;
  vector<vector<u64>> A(N + 1, vector<u64>(sz, 0));
  vector<vector<u64>> B(N + 1, vector<u64>(sz, 0));
  for (int i = 0; i < sz; ++i) {
   int p = __builtin_popcount(i);
   A[p][i] = a[i]; B[p][i] = b[i];
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
    for (int q = 0; q <= p; ++q) sum = (sum + A[q][i] * B[p - q][i]) % MOD;
    C[p][i] = sum;
   }
  }
  for (int k = 0; k <= N; ++k) mobius(C[k]);
  vector<u32> c(sz);
  for (int i = 0; i < sz; ++i) c[i] = (u32) C[__builtin_popcount(i)][i];
  return c;
 }
 static vector<u32> run(int N, int M, const vector<u32>& a, const vector<u32>& w) {
  int sz = 1 << N;
  vector<u32> b(sz, 0);
  b[0] = 1; // a^0 = identity
  vector<u32> ans(M, 0);
  for (int m = 0; m < M; ++m) {
   u64 s = 0;
   for (int i = 0; i < sz; ++i) s = (s + (u64) w[i] * b[i]) % MOD;
   ans[m] = (u32) s;
   if (m + 1 < M) b = conv_(N, b, a);
  }
  return ans;
 }
};
