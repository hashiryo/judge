#pragma once
// 自前ライブラリ (judge/lib/mylib) の sps::log<ModInt<998244353>> を呼ぶラッパー。
// O(n^2 · 2^n) の subset log。

#include "_common.hpp"
#include "mylib/algebra/ModInt.hpp"
#include "mylib/algebra/set_power_series.hpp"

struct SubsetLog {
 static vector<u32> run(int N, const vector<u32>& b_in) {
  using Mint = ModInt<998244353>;
  int sz = 1 << N;
  vector<Mint> f(sz);
  for (int i = 0; i < sz; ++i) f[i] = Mint(b_in[i]);
  auto c = sps::log(f);
  vector<u32> result(sz);
  for (int i = 0; i < sz; ++i) result[i] = (u32) c[i].val();
  return result;
 }
};
