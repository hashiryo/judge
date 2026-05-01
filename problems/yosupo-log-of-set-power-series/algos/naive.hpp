#pragma once
#include "_common.hpp"
// 半分再帰版 subset log: out0 = log(g[:N/2])
//                         out1 = g[N/2:] / g[:N/2]   (subset div)
//                         return concat(out0, out1)
// 中身の subset_div / subset_convolution は素朴 zeta/Möbius で O(M^2 · 2^M)。
struct SubsetLog {
 // 素朴 subset convolution
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
 static u32 qpow_(u64 x, u64 y) {
  u64 r = 1;
  while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
  return (u32) r;
 }
 // a / b: 共に nilpotent 部 (b[0] != 0 必須)。
 //   c[0] = a[0] / b[0]、それ以降は c = (a - c * (b - b[0])) / b[0]
 //   素朴に subset convolution + 引き算で漸進的に求める方が簡単 (ランクごと)。
 static vector<u32> div_(int N, const vector<u32>& a, const vector<u32>& b) {
  int sz = 1 << N;
  // rank vector で OR-convolution の zeta 領域へ。
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
  // C[k][i] = (A[k][i] - sum_{q=0..k-1} C[q][i] * B[k-q][i]) / B[0][i]
  vector<vector<u64>> C(N + 1, vector<u64>(sz, 0));
  for (int i = 0; i < sz; ++i) {
   u32 inv0 = qpow_(B[0][i], MOD - 2);
   for (int k = 0; k <= N; ++k) {
    u64 v = A[k][i];
    for (int q = 0; q < k; ++q) v = (v + MOD - C[q][i] * B[k - q][i] % MOD) % MOD;
    C[k][i] = v * inv0 % MOD;
   }
  }
  for (int k = 0; k <= N; ++k) mobius(C[k]);
  vector<u32> c(sz);
  for (int i = 0; i < sz; ++i) c[i] = (u32) C[__builtin_popcount(i)][i];
  return c;
 }
 static vector<u32> log_rec(int N, const vector<u32>& g) {
  if (N == 0) return {0}; // log(1) = 0
  int half = 1 << (N - 1);
  vector<u32> gl(g.begin(), g.begin() + half);
  vector<u32> gh(g.begin() + half, g.end());
  vector<u32> out0 = log_rec(N - 1, gl);
  vector<u32> out1 = div_(N - 1, gh, gl);
  vector<u32> out;
  out.reserve(1 << N);
  out.insert(out.end(), out0.begin(), out0.end());
  out.insert(out.end(), out1.begin(), out1.end());
  return out;
 }
 static vector<u32> run(int N, const vector<u32>& b_in) {
  return log_rec(N, b_in);
 }
};
