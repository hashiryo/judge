#pragma once
#include "_common.hpp"
// 教科書的 GF(2) ガウス消去。pivot 行数 = ランク。
// 制約: N*M ≤ 2^24 なので M は 2^24 まで取れる (N=1 のとき)。
// 行は u64 配列で動的サイズに保持 (bitset<4096> だと M>4096 で OOB → UB)。
struct Rank {
 static int run(int n, int m, const vector<string>& a) {
  const int W = (m + 63) / 64;
  vector<u64> M((size_t) n * W, 0);
  for (int i = 0; i < n; ++i) {
   const char* row = a[i].data();
   for (int j = 0; j < m; ++j)
    if (row[j] == '1') M[(size_t) i * W + j / 64] |= u64(1) << (j % 64);
  }
  int piv = 0;
  for (int j = 0; j < m && piv < n; ++j) {
   int sel = -1;
   for (int i = piv; i < n; ++i)
    if ((M[(size_t) i * W + j / 64] >> (j % 64)) & 1) { sel = i; break; }
   if (sel == -1) continue;
   if (sel != piv) {
    for (int w = 0; w < W; ++w) std::swap(M[(size_t) piv * W + w], M[(size_t) sel * W + w]);
   }
   for (int i = piv + 1; i < n; ++i) {
    if ((M[(size_t) i * W + j / 64] >> (j % 64)) & 1) {
     for (int w = 0; w < W; ++w) M[(size_t) i * W + w] ^= M[(size_t) piv * W + w];
    }
   }
   ++piv;
  }
  return piv;
 }
};
