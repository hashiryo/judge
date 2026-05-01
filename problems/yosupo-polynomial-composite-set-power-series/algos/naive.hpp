#pragma once
#include "_common.hpp"
// 素朴: f(b) = sum_{j=0..M-1} f[j] * b^j。
//   acc = b^0 = identity (acc[0]=1)
//   res += f[0] * acc; for j=1..M-1: acc = acc * b; res += f[j] * acc
//   subset_convolution は素朴 zeta/Möbius で O(N^2 · 2^N)。
//   全体で O(M · N^2 · 2^N)。M=1e5, N=20 だと TLE 必至。
struct PolyComp {
 static vector<u32> conv_(int N, const vector<u32>& a, const vector<u32>& b) {
  int sz = 1 << N;
  vector<vector<u64>> A(N + 1, vector<u64>(sz, 0));
  vector<vector<u64>> B(N + 1, vector<u64>(sz, 0));
  for (int i = 0; i < sz; ++i) {
   int p = __builtin_popcount(i);
   A[p][i] = a[i]; B[p][i] = b[i];
  }
  auto zeta = [&](vector<u64>& F) {
   for (int i = 0; i < N; ++i)
    for (int j = 0; j < sz; ++j)
     if (j & (1 << i)) F[j] = (F[j] + F[j ^ (1 << i)]) % MOD;
  };
  auto mobius = [&](vector<u64>& F) {
   for (int i = 0; i < N; ++i)
    for (int j = 0; j < sz; ++j)
     if (j & (1 << i)) F[j] = (F[j] + MOD - F[j ^ (1 << i)]) % MOD;
  };
  for (int k = 0; k <= N; ++k) zeta(A[k]);
  for (int k = 0; k <= N; ++k) zeta(B[k]);
  vector<vector<u64>> C(N + 1, vector<u64>(sz, 0));
  for (int i = 0; i < sz; ++i) {
   for (int p = 0; p <= N; ++p) {
    u64 s = 0;
    for (int q = 0; q <= p; ++q) s = (s + A[q][i] * B[p - q][i]) % MOD;
    C[p][i] = s;
   }
  }
  for (int k = 0; k <= N; ++k) mobius(C[k]);
  vector<u32> c(sz);
  for (int i = 0; i < sz; ++i) c[i] = (u32) C[__builtin_popcount(i)][i];
  return c;
 }
 static vector<u32> run(int M, int N, const vector<u32>& f, const vector<u32>& b) {
  int sz = 1 << N;
  vector<u32> res(sz, 0);
  if (M == 0) return res;
  vector<u32> acc(sz, 0); acc[0] = 1; // b^0 = identity
  for (int j = 0; j < M; ++j) {
   if (f[j]) {
    for (int i = 0; i < sz; ++i)
     res[i] = (u32) ((res[i] + (u64) f[j] * acc[i]) % MOD);
   }
   if (j + 1 < M) acc = conv_(N, acc, b);
  }
  return res;
 }
};
