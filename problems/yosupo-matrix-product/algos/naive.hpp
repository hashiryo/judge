#pragma once
#include "_common.hpp"
// O(N*M*P) 素朴行列積。mod 演算は最後にまとめて。
struct Mul {
 static vector<u32> run(int n, int m, int p, const vector<u32>& a, const vector<u32>& b) {
  vector<u32> c((size_t) n * p, 0);
  for (int i = 0; i < n; ++i) {
   for (int k = 0; k < m; ++k) {
    u64 aik = a[(size_t) i * m + k];
    if (!aik) continue;
    for (int j = 0; j < p; ++j) {
     c[(size_t) i * p + j] = (u32) ((c[(size_t) i * p + j] + aik * b[(size_t) k * p + j]) % MOD);
    }
   }
  }
  return c;
 }
};
