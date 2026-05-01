#pragma once
// =============================================================================
// Source: yosupo "Matrix Rank (mod 2)" 提出 345398 を抽出。
//   https://judge.yosupo.jp/submission/345398
//   方式: 普通のガウス消去 + AVX2 で XOR を 4 命令ずつアンロール。
//   N < M なら入力を転置して扱う (rank 不変、cols が小さい方が速い)。
// 抽出方針:
//   - I/O は base.cpp 経由
//   - Matrix struct はそのまま、AVX2 経路は USE_SIMDE 環境でも有効化
//   - Rank::run wrapper で string ↔ Matrix セット
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

// !!! _common.hpp より先に simde を include (float16_t 衝突回避)
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

namespace yosupo_345398 {
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
 bool get(int r, int c) { return (data[r * stride + (c / 64)] >> (c % 64)) & 1; }
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
} // namespace yosupo_345398

struct Rank {
 static int run(int n, int m, const vector<string>& a) {
  using namespace yosupo_345398;
  // N < M なら転置すると cols が小さくなり高速化。rank は不変。
  if (n < m) {
   Matrix mat(n, m);
   for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
     if (a[i][j] == '1') mat.set(i, j);
   return mat.gaussian_elimination();
  } else {
   Matrix mat(m, n);
   for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
     if (a[i][j] == '1') mat.set(j, i);
   return mat.gaussian_elimination();
  }
 }
};
