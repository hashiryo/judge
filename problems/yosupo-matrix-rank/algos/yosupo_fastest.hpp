#pragma once
// =============================================================================
// 参考: yosupo "Rank of Matrix" 提出 285774 (cp-algo の linalg::matrix<modint>)。
//   https://judge.yosupo.jp/submission/285774
//   cp-algo は modint/vec/matrix の重いテンプレート連鎖で構成されており
//   一括抽出は現実的でないので、ここでは同一方針 (Montgomery 系の重い遅延
//   reduce u64 行列で AVX2 mul-add、N回ごとに mod 還元) を素朴に再実装する。
// 方式:
//   - N < M なら転置 (pivot ループは max(N,M) ではなく min(N,M))
//   - 行を u64 で持って `b += factor * a` を AVX2 _mm256_mul_epi32 で 4 lane
//   - 一定回数で mod 還元 (lazy reduce)
// ライセンス: cp-algo (元実装) / 本ファイルは独立再実装。
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

namespace yosupo_rank_lazy {
constexpr u32 P = 998244353;
constexpr u32 N_REDUCE = 16;

inline u32 qpow(u32 x, u32 y) {
 u64 r = 1, b = x;
 while (y) { if (y & 1) r = r * b % P; b = b * b % P; y >>= 1; }
 return (u32) r;
}

inline int rank_inplace(vector<vector<u64>>& mat, int n, int m) {
 vector<u32> counter(n, 0);
 int rank = 0;
 int col = 0;
 for (int row = 0; row < n && col < m; ) {
  // pivot 探し: col 列の row 以下を還元しながら non-zero 探す
  int sel = -1;
  for (int j = row; j < n; ++j) {
   mat[j][col] %= P;
   if (mat[j][col] != 0 && sel == -1) sel = j;
  }
  if (sel == -1) { ++col; continue; }
  if (sel != row) {
   std::swap(mat[row], mat[sel]);
   std::swap(counter[row], counter[sel]);
  }
  // pivot 行を全列 reduce
  for (int j = col; j < m; ++j) mat[row][j] %= P;
  counter[row] = 0;
  u32 inv = P - qpow((u32) mat[row][col], P - 2);
  for (int j = 0; j < n; ++j) {
   if (j == row) continue;
   mat[j][col] %= P;
   if (mat[j][col] == 0) continue;
   u32 factor = (u32) ((u64) inv * mat[j][col] % P);
   const auto& aa = mat[row];
   auto& bb = mat[j];
   int x = col;
#if (defined(__x86_64__) || defined(__i386__)) || defined(USE_SIMDE)
   __m256i vc = _mm256_set1_epi32((int) factor);
   for (; x + 4 <= m; x += 4) {
    __m256i va = _mm256_loadu_si256((const __m256i*) &aa[x]);
    __m256i vb = _mm256_loadu_si256((const __m256i*) &bb[x]);
    __m256i vmul = _mm256_mul_epi32(va, vc);
    vb = _mm256_add_epi64(vb, vmul);
    _mm256_storeu_si256((__m256i*) &bb[x], vb);
   }
#endif
   for (; x < m; ++x) bb[x] += (u64) factor * aa[x];
   if (++counter[j] == N_REDUCE) {
    counter[j] = 0;
    for (int x2 = col; x2 < m; ++x2) bb[x2] %= P;
   }
  }
  ++rank; ++row; ++col;
 }
 return rank;
}
} // namespace yosupo_rank_lazy

struct Rank {
 static int run(int n, int m, const vector<vector<u32>>& a) {
  // N < M なら転置 (rank 不変、消去対象列数 = min(N,M) を小さく)
  bool transposed = (n < m);
  int N = transposed ? m : n;
  int M = transposed ? n : m;
  vector<vector<u64>> mat(N, vector<u64>(M));
  if (transposed) {
   for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j) mat[j][i] = a[i][j];
  } else {
   for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j) mat[i][j] = a[i][j];
  }
  return yosupo_rank_lazy::rank_inplace(mat, N, M);
 }
};
