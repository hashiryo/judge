#pragma once
#include "_common.hpp"
// =============================================================================
// 比較用: 任意 mod 行列式の fastest を 998244353 固定で動かしたもの。
//   元提出 (matrix_det_arbitrary_mod #361230, 2 位): https://judge.yosupo.jp/submission/361230
//   方式: Barrett64 reduction + 互除法 Gauss。素数性を仮定しないので素数特化版より遅いはず。
//
// 比較ポイント:
//   - 同じ行列式問題でも「素数法と分かっていれば逆元を qpow で出せる」のと
//     「任意 mod なので互除法 Gauss」では計算量が違う (互除法は log M 倍重い)
//   - 素数 998244353 でこの algo を動かすと「不要な汎用性のオーバヘッド」が出る
//   - 本来の det 高速版 (yosupo_fastest.hpp = #361510) が prime branch + AVX2 で
//     どれだけ得をしているかが定量化できる
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")

namespace yosupo_361230_for_det {
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
    if (j != i) { std::swap(mat[i], mat[j]); neg ^= 1; }
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
} // namespace yosupo_361230_for_det

struct Det {
 static u32 run(int n, const vector<vector<u32>>& a) {
  vector<vector<u32>> mat = a;
  return yosupo_361230_for_det::det_inplace(mat, (u32) n, MOD);
 }
};
