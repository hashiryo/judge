#pragma once
// =============================================================================
// Source: yosupo "Inverse Matrix (mod 2)" 提出 241305 を抽出。
//   https://judge.yosupo.jp/submission/241305
//   方式: [A | I] (N×2N) を作り、16-row block のブロック化ガウス・ジョルダン消去。
//   - 前進消去: vecs[i] の pivot 列を firsts[i] に記録 (行 swap せず)
//   - 後退消去: 上三角部分を消して I を逆行列に変換
//   - branchless XOR: vecs[j] と zeroes の三項選択で常に xor_from
//   出力時に invfirsts[firsts[i]] = i で行順を復元。
// 抽出方針:
//   - I/O は base.cpp 経由
//   - グローバル配列は namespace yosupo_241305 に閉じる
//   - Inv::run wrapper で string ↔ ビットパック変換、出力で右半分を抽出
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

// !!! _common.hpp より先に simde を include (float16_t 衝突回避)
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

namespace yosupo_241305 {
constexpr int FN = 4096;
constexpr int LEN = FN / 256 * 2;  // = 32 (左 N + 右 N の半分ずつ)
constexpr int BLOCKSIZE = 16;

inline uint16_t firsts[FN];
inline uint16_t invfirsts[FN];
inline __m256i zeroes[LEN];
inline __m256i vecs[FN][LEN];

inline void xor_from(__m256i* dst, const __m256i* src, size_t start) {
 for (size_t i = start / 256; i < LEN; ++i) dst[i] ^= src[i];
}
inline size_t find_first_before(__m256i* vec, size_t n) {
 uint64_t* v64 = (uint64_t*) vec;
 for (size_t i = 0; i < (n + 63) / 64; ++i) {
  if (v64[i] != 0) {
   size_t ans = i * 64 + (size_t) __builtin_ctzll(v64[i]);
   return ans < n ? ans : ~(size_t) 0;
  }
 }
 return ~(size_t) 0;
}
inline int get_ith(const __m256i* vec, size_t idx) {
 return (((const uint64_t*) vec)[idx / 64] >> (idx % 64)) & 1;
}
inline void set_ith(__m256i* vec, size_t idx) {
 ((uint64_t*) vec)[idx / 64] |= uint64_t(1) << (idx % 64);
}
} // namespace yosupo_241305

struct Inv {
 static vector<string> run(int n_in, const vector<string>& a) {
  using namespace yosupo_241305;
  std::memset(firsts, 0, sizeof(firsts));
  std::memset(invfirsts, 0, sizeof(invfirsts));
  std::memset(zeroes, 0, sizeof(zeroes));
  std::memset(vecs, 0, sizeof(vecs));
  // a を vecs[i] の左半分にパック。
  for (int i = 0; i < n_in; ++i) {
   const char* row = a[i].data();
   uint32_t* d32 = (uint32_t*) vecs[i];
   for (int j = 0; j < n_in; ++j)
    if (row[j] == '1') d32[j / 32] |= uint32_t(1) << (j % 32);
  }
  size_t orig_n = (size_t) n_in;
  size_t n = orig_n;
  // n を BLOCKSIZE の倍数に揃える (恒等行を追加するだけで det 不変)。
  for (; n % BLOCKSIZE != 0; ++n) set_ith(vecs[n], n);
  // 右半分に単位行列をセット (列 n..n+n-1)。
  for (size_t i = 0; i < n; ++i) set_ith(vecs[i], i + n);

  // 前進消去
  for (size_t il = 0; il < n; il += BLOCKSIZE) {
   for (size_t jl = 0; jl < il; jl += BLOCKSIZE) {
    for (size_t j = jl; j < jl + BLOCKSIZE; ++j)
     for (size_t i = il; i < il + BLOCKSIZE; ++i) {
      __m256i* src = get_ith(vecs[i], firsts[j]) ? vecs[j] : zeroes;
      xor_from(vecs[i], src, firsts[j]);
     }
   }
   for (size_t i = il; i < il + BLOCKSIZE; ++i) {
    for (size_t j = il; j < i; ++j) {
     __m256i* src = get_ith(vecs[i], firsts[j]) ? vecs[j] : zeroes;
     xor_from(vecs[i], src, firsts[j]);
    }
    size_t first_bit = find_first_before(vecs[i], n);
    if (first_bit == ~(size_t) 0) return {};  // 可逆でない
    firsts[i] = (uint16_t) first_bit;
   }
  }

  // 後退消去
  for (size_t jl = 0; jl < n; jl += BLOCKSIZE) {
   for (size_t il = jl + BLOCKSIZE; il < n; il += BLOCKSIZE) {
    for (size_t i = il; i < il + BLOCKSIZE; ++i)
     for (size_t j = jl; j < jl + BLOCKSIZE; ++j) {
      __m256i* src = get_ith(vecs[j], firsts[i]) ? vecs[i] : zeroes;
      xor_from(vecs[j], src, firsts[i]);
     }
   }
   for (size_t j = jl; j < jl + BLOCKSIZE; ++j) {
    for (size_t i = j + 1; i < jl + BLOCKSIZE; ++i) {
     __m256i* src = get_ith(vecs[j], firsts[i]) ? vecs[i] : zeroes;
     xor_from(vecs[j], src, firsts[i]);
    }
   }
  }

  // 行順を invfirsts で復元しつつ、右半分 [n, n + orig_n) を出力に変換。
  for (size_t i = 0; i < n; ++i) invfirsts[firsts[i]] = (uint16_t) i;
  vector<string> out(orig_n, string(orig_n, '0'));
  for (size_t i = 0; i < orig_n; ++i) {
   const __m256i* row = vecs[invfirsts[i]];
   const uint64_t* r64 = (const uint64_t*) row;
   for (size_t j = 0; j < orig_n; ++j) {
    size_t bit = j + n;
    if ((r64[bit / 64] >> (bit % 64)) & 1) out[i][j] = '1';
   }
  }
  return out;
 }
};
