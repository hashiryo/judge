#pragma once
// =============================================================================
// Source: yosupo "Determinant of Matrix" 提出 361510 を抽出。
//   https://judge.yosupo.jp/submission/361510
//   方式: u64 行列で「N=16 回の積和に 1 回だけ mod 還元」する遅延 reduce。
//   AVX2 `_mm256_mul_epi32` で 4 lane × 64-bit を一気に積算 (lane あたり 32-bit
//   入力 → 64-bit 結果)。Counter で各行の累積回数を追跡、N に達したら mod 取る。
// 抽出方針:
//   - I/O 部分 (FastIO) は base.cpp に移譲して削除
//   - prime branch のみ抽出 (mod 998244353 用)
//   - Det::run wrapper で vector<vector<u32>> → vector<vector<u64>> 変換
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

namespace yosupo_361510 {
constexpr u32 P = 998244353;
constexpr u32 N_REDUCE = 16;  // 16 回積和ごとに 1 回 mod 還元

inline u32 qpow(u32 x, u32 y) {
 u64 res = 1, b = x;
 while (y) { if (y & 1) res = res * b % P; b = b * b % P; y >>= 1; }
 return (u32) res;
}

inline u32 det_inplace(vector<vector<u64>>& mat, u32 n) {
 if (n == 0) return 1;
 bool neg = false;
 vector<u32> counter(n, 0);
 for (u32 i = 0; i < n; ++i) {
  for (u32 j = i; j < n; ++j) mat[j][i] %= P;
  for (u32 j = i; j < n; ++j) {
   if (mat[j][i] != 0) {
    if (j != i) {
     std::swap(mat[i], mat[j]);
     counter[j] = counter[i];
     neg ^= 1;
    }
    break;
   }
  }
  for (u32 j = i; j < n; ++j) mat[i][j] %= P;
  if (mat[i][i] == 0) return 0;
  u32 inv = P - qpow((u32) mat[i][i], P - 2);
  for (u32 j = i + 1; j < n; ++j) {
   u32 factor = (u32) ((u64) inv * mat[j][i] % P);
   const auto& a = mat[i];
   auto& b = mat[j];
   u32 x = i;
   __m256i vc = _mm256_set1_epi32((int) factor);
   for (; x + 4 <= n; x += 4) {
    __m256i va = _mm256_loadu_si256((const __m256i*) &a[x]);
    __m256i vb = _mm256_loadu_si256((const __m256i*) &b[x]);
    __m256i vmul = _mm256_mul_epi32(va, vc);
    vb = _mm256_add_epi64(vb, vmul);
    _mm256_storeu_si256((__m256i*) &b[x], vb);
   }
   for (; x < n; ++x) b[x] += (u64) factor * a[x];
   if (++counter[j] == N_REDUCE) {
    counter[j] = 0;
    for (u32 x2 = i; x2 < n; ++x2) b[x2] %= P;
   }
  }
 }
 u64 res = neg ? P - 1 : 1;
 for (u32 i = 0; i < n; ++i) res = res * mat[i][i] % P;
 return (u32) res;
}
} // namespace yosupo_361510

struct Det {
 static u32 run(int n, const vector<vector<u32>>& a) {
  vector<vector<u64>> mat(n, vector<u64>(n));
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j) mat[i][j] = a[i][j];
  return yosupo_361510::det_inplace(mat, (u32) n);
 }
};
