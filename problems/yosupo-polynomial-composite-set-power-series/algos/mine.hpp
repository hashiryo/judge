#pragma once
// 自前ライブラリ (judge/lib/mylib) の sps::poly_comp<ModInt<998244353>> を呼ぶ。
// O(n^2 · 2^n + ?)。

#include "_common.hpp"
#include "mylib/algebra/ModInt.hpp"
#include "mylib/algebra/set_power_series.hpp"

struct PolyComp {
 static vector<u32> run(int M, int N, const vector<u32>& f_in, const vector<u32>& b_in) {
  using Mint = ModInt<998244353>;
  int sz = 1 << N;
  vector<Mint> P(M);
  for (int i = 0; i < M; ++i) P[i] = Mint(f_in[i]);
  vector<Mint> ff(sz);
  for (int i = 0; i < sz; ++i) ff[i] = Mint(b_in[i]);
  auto c = sps::poly_comp(P, ff);
  vector<u32> result(sz);
  for (int i = 0; i < sz; ++i) result[i] = (u32) c[i].val();
  return result;
 }
};
