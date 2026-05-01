#pragma once
// =============================================================================
// 自前ライブラリ (judge/lib/mylib = algo/Library サブモジュール) の
//   sps::convolve<ModInt<998244353>> を呼ぶラッパー。
// 中身は rank zeta / Möbius (nibble-packed rank-vector + zeta) を用いた
// O(N^2 2^N) subset convolution。N < 11 の小さいケースは naive に切り替わる。
// 比較用なので最速の SIMD は意識していない。
// =============================================================================

#include "_common.hpp"
#include "mylib/algebra/ModInt.hpp"
#include "mylib/algebra/set_power_series.hpp"

struct SubsetConv {
 static vector<u32> run(int N, const vector<u32>& a_in, const vector<u32>& b_in) {
  using Mint = ModInt<998244353>;
  int sz = 1 << N;
  vector<Mint> a(sz), b(sz);
  for (int i = 0; i < sz; ++i) {
   a[i] = Mint(a_in[i]);
   b[i] = Mint(b_in[i]);
  }
  auto c = sps::convolve(a, b);
  vector<u32> result(sz);
  for (int i = 0; i < sz; ++i) result[i] = (u32) c[i].val();
  return result;
 }
};
