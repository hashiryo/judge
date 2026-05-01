#pragma once
#include "_common.hpp"
// =============================================================================
// Source: yosupo "Matrix Determinant (arbitrary mod)" 提出 361230 を抽出 (2 位)。
//   https://judge.yosupo.jp/submission/361230
//   方式: Barrett reduction (u64) で mod 演算を高速化した素朴 Gauss 消去。
//   AVX-512 不要、portable。N ≤ 500 で十分速い。
// 抽出方針:
//   - I/O は base.cpp 経由 (元コードの fast I/O は不要)
//   - Barrett64 + det_inplace を namespace に閉じる
//   - Det::run wrapper で N, mod, A → u32 結果
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")

namespace yosupo_361230 {
struct Barrett64 {
 uint64_t mod, inv;
 Barrett64() = default;
 explicit Barrett64(uint64_t m): mod(m), inv(-1ull / m + 1) {}
 constexpr uint64_t reduce(uint64_t x) const {
  x -= static_cast<uint64_t>(static_cast<__uint128_t>(x) * inv >> 64) * mod;
  x += -(x >> 63) & mod;
  return x;
 }
};

inline uint32_t det_inplace(vector<vector<u32>>& mat, uint32_t n, uint32_t P) {
 if (P == 1) return 0;
 Barrett64 bar(P);
 bool neg = false;
 for (uint32_t i = 0; i < n; ++i) {
  for (uint32_t j = i; j < n; ++j) {
   if (mat[j][i] != 0) {
    if (j != i) {
     std::swap(mat[i], mat[j]);
     neg ^= 1;
    }
    break;
   }
  }
  for (uint32_t j = i + 1; j < n; ++j) {
   while (mat[i][i] != 0) {
    uint32_t inv = P - mat[j][i] / mat[i][i];
    if (inv != P) {
     for (uint32_t x = i; x < n; ++x) {
      mat[j][x] = (uint32_t) bar.reduce(mat[j][x] + 1ull * inv * mat[i][x]);
     }
    }
    std::swap(mat[i], mat[j]);
    neg ^= 1;
   }
   std::swap(mat[i], mat[j]);
   neg ^= 1;
  }
 }
 uint32_t res = neg ? P - 1 : 1;
 for (uint32_t i = 0; i < n; ++i) res = (uint32_t) bar.reduce((uint64_t) res * mat[i][i]);
 return res;
}
} // namespace yosupo_361230

struct Det {
 static u32 run(int n, u32 mod, const vector<vector<u32>>& a) {
  vector<vector<u32>> mat = a;
  return yosupo_361230::det_inplace(mat, (u32) n, mod);
 }
};
