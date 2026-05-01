#pragma once
#include "_common.hpp"
// 標準的な O(N^2 2^N) subset exp:
//   exp は subset_convolution の漸化式で計算できる。
//     out[empty] = 1
//     for k in 1..N:
//       out[low | {k}] = sum over s in low: b[s | {k}] * out[low \ s] * (...)/k
//   最も簡単な実装: 半分ずつ再帰
//     N==0 → {1}
//     out0 = subset_exp(b の lower half (size 2^(N-1)))
//     out1 = subset_convolution(out0, b の upper half)
//     return concat(out0, out1)
//   ただし subset_convolution は素朴 zeta/Möbius で O(M^2 · 2^M)。
//   全体で T(N) = T(N-1) + O(N^2 · 2^N) → O(N^2 · 2^N)。
struct SubsetExp {
 // 素朴 subset convolution (rank vector + zeta/Möbius)
 static vector<u32> conv_(int N, const vector<u32>& a, const vector<u32>& b) {
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
    for (int q = 0; q <= p; ++q) sum = (sum + A[q][i] * B[p - q][i]) % MOD;
    C[p][i] = sum;
   }
  }
  for (int k = 0; k <= N; ++k) mobius(C[k]);
  vector<u32> c(sz);
  for (int i = 0; i < sz; ++i) c[i] = (u32) C[__builtin_popcount(i)][i];
  return c;
 }
 static vector<u32> exp_rec(int N, const vector<u32>& b) {
  if (N == 0) return {1};
  int half = 1 << (N - 1);
  vector<u32> bl(b.begin(), b.begin() + half);
  vector<u32> bh(b.begin() + half, b.end());
  vector<u32> out0 = exp_rec(N - 1, bl);
  vector<u32> out1 = conv_(N - 1, out0, bh);
  vector<u32> out;
  out.reserve(1 << N);
  out.insert(out.end(), out0.begin(), out0.end());
  out.insert(out.end(), out1.begin(), out1.end());
  return out;
 }
 static vector<u32> run(int N, const vector<u32>& b_in) {
  return exp_rec(N, b_in);
 }
};
