#pragma once
// 自前ライブラリ (judge/lib/mylib) の is_prime を呼ぶラッパー。
// 範囲別に Montgomery 32 / Montgomery 64 / Barrett 系を切り替えつつ
// Miller-Rabin (固定 base) で判定する。

#include "_common.hpp"
#include "mylib/number_theory/is_prime.hpp"

struct Primality {
 static vector<bool> run(const vector<u64>& qs) {
  vector<bool> ans(qs.size());
  for (size_t i = 0; i < qs.size(); ++i) ans[i] = is_prime(qs[i]);
  return ans;
 }
};
