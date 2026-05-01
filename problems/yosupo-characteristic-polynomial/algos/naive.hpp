#pragma once
#include "_common.hpp"
// Faddeev-LeVerrier 法 (O(N^4)):
//   M_0 = I, c_N = 1
//   for k = 1..N:
//     M_k = M * (M_{k-1} - c_{N-k+1} I)  (これは k 回目の中間行列)
//     c_{N-k} = -tr(M * M_{k-1}) / k
//   結果として p(x) = sum c_i x^i, c_N = 1。
// 998244353 は素数で N ≤ 500 < MOD なので k は逆元あり。
struct CharPoly {
 static u32 qpow_(u64 x, u64 y) {
  u64 r = 1;
  while (y) { if (y & 1) r = r * x % MOD; x = x * x % MOD; y >>= 1; }
  return (u32) r;
 }
 static vector<u32> run(int n, const vector<vector<u32>>& M_in) {
  vector<u32> c(n + 1, 0);
  c[n] = 1;
  if (n == 0) return c;
  // Mk: 現在の M^k 様の行列。最初は M。一般に M_{k} (上の漸化式) を保持。
  vector<vector<u32>> Mk(n, vector<u32>(n, 0));
  for (int i = 0; i < n; ++i) Mk[i][i] = 1; // M_0 = I
  // 1 step: Mk_new = M * (Mk - c[n-k+1] I)、ただし Mk は M_{k-1}
  vector<vector<u32>> tmp(n, vector<u32>(n, 0));
  for (int k = 1; k <= n; ++k) {
   // tmp = M * Mk
   for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
     u64 s = 0;
     for (int l = 0; l < n; ++l)
      s += (u64) M_in[i][l] * Mk[l][j] % MOD;
     tmp[i][j] = (u32) (s % MOD);
    }
   }
   // c[n-k] = -tr(tmp) / k
   u64 tr = 0;
   for (int i = 0; i < n; ++i) tr = (tr + tmp[i][i]) % MOD;
   u32 inv_k = qpow_(k, MOD - 2);
   u32 ck = (u32) ((u64) (MOD - tr) % MOD * inv_k % MOD);
   c[n - k] = ck;
   // Mk_new = tmp + c[n-k] * I  (= M * (M_{k-1} - c_{N-k+1} I) ではないが、
   //  Faddeev-LeVerrier の正しい漸化式は M_k = M * M_{k-1} + c_{N-k} I)
   if (k < n) {
    for (int i = 0; i < n; ++i)
     for (int j = 0; j < n; ++j)
      Mk[i][j] = tmp[i][j];
    for (int i = 0; i < n; ++i)
     Mk[i][i] = (u32) ((Mk[i][i] + ck) % MOD);
   }
  }
  return c;
 }
};
