#pragma once
// =============================================================================
// 「rank の計算アルゴリズムで det を解く」変種。
//   GF(2) では det(A) = 1 ⟺ rank(A) = N なので、rank の問題の最速解 (#345398)
//   をそのまま使って rank == N を判定するだけで det が出る。
//
//   元コード: yosupo "Matrix Rank (mod 2)" 提出 345398
//     https://judge.yosupo.jp/submission/345398
//
//   比較ポイント:
//   - dedicated 版 (#241296) は 32-row block の branchless Gauss
//   - rank 経由 (#345398) は素直なガウス + AVX2 XOR 4-アンロール
//   det では行列が常に N×N なので、rank 版の「N<M なら転置」トリックは不要。
//   非可逆 (det=0) の場合、dedicated 版は pivot 失敗時に即 return できるが、
//   rank 版は最後まで走る (rank が N 未満かどうかを最後まで決定するため)
//   ── ただし元の gaussian_elimination は「pivot 見つからなければ列スキップ」
//   して continue するだけなので、無駄に時間を使う訳ではない。
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#define USE_AVX2
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#define USE_AVX2
#endif

#include "_common.hpp"

namespace yosupo_345398_det_via_rank {
struct Matrix {
 const int rows, cols, stride;
 unsigned long long* data;
 Matrix(int r, int c): rows(r), cols(c), stride((c + 63) / 64) {
  size_t total_bytes = (size_t) rows * stride * sizeof(unsigned long long);
  data = new unsigned long long[(size_t) rows * stride];
  std::memset(data, 0, total_bytes);
 }
 ~Matrix() { delete[] data; }
 void set(int r, int c) { data[r * stride + (c / 64)] |= (1ULL << (c % 64)); }
#ifdef USE_AVX2
 void xor_rows(unsigned long long* d, const unsigned long long* s, int start_idx) {
  int i = start_idx;
  for (; i + 16 <= stride; i += 16) {
   _mm256_storeu_si256((__m256i*) &d[i], _mm256_xor_si256(_mm256_loadu_si256((__m256i*) &d[i]), _mm256_loadu_si256((__m256i*) &s[i])));
   _mm256_storeu_si256((__m256i*) &d[i + 4], _mm256_xor_si256(_mm256_loadu_si256((__m256i*) &d[i + 4]), _mm256_loadu_si256((__m256i*) &s[i + 4])));
   _mm256_storeu_si256((__m256i*) &d[i + 8], _mm256_xor_si256(_mm256_loadu_si256((__m256i*) &d[i + 8]), _mm256_loadu_si256((__m256i*) &s[i + 8])));
   _mm256_storeu_si256((__m256i*) &d[i + 12], _mm256_xor_si256(_mm256_loadu_si256((__m256i*) &d[i + 12]), _mm256_loadu_si256((__m256i*) &s[i + 12])));
  }
  for (; i < stride; ++i) d[i] ^= s[i];
 }
#else
 void xor_rows(unsigned long long* d, const unsigned long long* s, int start_idx) {
  for (int i = start_idx; i < stride; ++i) d[i] ^= s[i];
 }
#endif
 int gaussian_elimination() {
  int pivot_row = 0;
  for (int col = 0; col < cols && pivot_row < rows; ++col) {
   int chunk_idx = col / 64;
   unsigned long long bit_mask = (1ULL << (col % 64));
   int sel = -1;
   for (int i = pivot_row; i < rows; ++i) {
    if (data[i * stride + chunk_idx] & bit_mask) { sel = i; break; }
   }
   if (sel == -1) continue;
   if (pivot_row != sel) xor_rows(data + pivot_row * stride, data + sel * stride, chunk_idx);
   for (int i = pivot_row + 1; i < rows; ++i) {
    if (data[i * stride + chunk_idx] & bit_mask) xor_rows(data + i * stride, data + pivot_row * stride, chunk_idx);
   }
   ++pivot_row;
  }
  return pivot_row;
 }
};
} // namespace yosupo_345398_det_via_rank

struct Det {
 static int run(int n, const vector<string>& a) {
  using namespace yosupo_345398_det_via_rank;
  Matrix mat(n, n);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j)
    if (a[i][j] == '1') mat.set(i, j);
  // rank == n ⟺ det = 1
  return mat.gaussian_elimination() == n ? 1 : 0;
 }
};
