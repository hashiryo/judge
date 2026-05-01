#pragma once
// =============================================================================
// Source: yosupo "Matrix Determinant (mod 2)" 提出 241296 を抽出。
//   https://judge.yosupo.jp/submission/241296
//   方式: GF(2) ガウス消去を 32-row block で再構成。pivot 行を swap せず
//   firsts[j] で各行の pivot 列を覚え、zeroes ベクトルとの三項演算
//   (vecs[j] : zeroes) で分岐を消した branchless XOR。
//   N=4096 固定で行列を 256-bit 単位 (NWORDS=16) でパック。
// 抽出方針:
//   - I/O は base.cpp 経由
//   - グローバル配列は namespace yosupo_241296 に閉じる
//   - Det::run wrapper で string ↔ ビットパック変換
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

// !!! _common.hpp より先に simde を include する !!!
// using namespace std で std::float16_t がグローバルに持ち込まれると aarch64 の
// arm_neon.h と衝突するため。
#pragma GCC optimize("O3,unroll-loops,rename-registers")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2,bmi")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#include <simde/x86/bmi.h>
#else
#include <immintrin.h>
#endif

#include "_common.hpp"

namespace yosupo_241296 {
constexpr int FN = 4096;
constexpr int NWORDS = FN / 256;  // = 16
constexpr int BLOCKSIZE = 32;

inline uint16_t firsts[FN];
inline __m256i zeroes[NWORDS];
inline __m256i vecs[FN][NWORDS];

inline void xor_from(__m256i* dst, const __m256i* src, size_t start) {
 for (size_t i = start / 256; i < NWORDS; ++i) dst[i] ^= src[i];
}
inline size_t find_first(__m256i* vec) {
 uint64_t* v64 = (uint64_t*) vec;
 for (size_t i = 0; i < FN / 64; ++i) {
  if (v64[i] != 0) return i * 64 + (size_t) __builtin_ctzll(v64[i]);
 }
 return ~(size_t) 0;
}
inline int get_ith(const __m256i* vec, size_t idx) {
 return (((const uint64_t*) vec)[idx / 64] >> (idx % 64)) & 1;
}
inline void set_ith(__m256i* vec, size_t idx) {
 ((uint64_t*) vec)[idx / 64] |= uint64_t(1) << (idx % 64);
}
} // namespace yosupo_241296

struct Det {
 static int run(int n_in, const vector<string>& a) {
  using namespace yosupo_241296;
  // 毎呼び出しで初期化。
  std::memset(firsts, 0, sizeof(firsts));
  std::memset(zeroes, 0, sizeof(zeroes));
  std::memset(vecs, 0, sizeof(vecs));
  // a[i][j] (char '0'/'1') を vecs[i] のビットにパック。
  for (int i = 0; i < n_in; ++i) {
   const char* row = a[i].data();
   uint32_t* d32 = (uint32_t*) vecs[i];
   for (int j = 0; j < n_in; ++j)
    if (row[j] == '1') d32[j / 32] |= uint32_t(1) << (j % 32);
  }
  size_t n = (size_t) n_in;
  // n を BLOCKSIZE の倍数に揃える (恒等行を追加するだけなので det は変わらない)。
  for (; n % BLOCKSIZE != 0; ++n) set_ith(vecs[n], n);

  for (size_t il = 0; il < n; il += BLOCKSIZE) {
   for (size_t jl = 0; jl < il; jl += BLOCKSIZE) {
    for (size_t j = jl; j < jl + BLOCKSIZE; ++j) {
     for (size_t i = il; i < il + BLOCKSIZE; ++i) {
      __m256i* src = get_ith(vecs[i], firsts[j]) ? vecs[j] : zeroes;
      xor_from(vecs[i], src, firsts[j]);
     }
    }
   }
   for (size_t i = il; i < il + BLOCKSIZE; ++i) {
    for (size_t j = il; j < i; ++j) {
     __m256i* src = get_ith(vecs[i], firsts[j]) ? vecs[j] : zeroes;
     xor_from(vecs[i], src, firsts[j]);
    }
    size_t first_bit = find_first(vecs[i]);
    if (first_bit == ~(size_t) 0) return 0;
    firsts[i] = (uint16_t) first_bit;
   }
  }
  return 1;
 }
};
