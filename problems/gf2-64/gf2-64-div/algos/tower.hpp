#pragma once
// 塔分解アプローチ: poly → nim → Nimber 演算 → poly。
// div: a / b = a * b^{-1}。Nimber::operator/ または Nimber::inv を使う。
#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/basis_change.hpp"
#include "mylib/algebra/Nimber.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_basis::poly_to_nim;
  using gf2_64_basis::nim_to_poly;
  Nimber::init();
  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   u64 a_nim = poly_to_nim(as[i]);
   u64 b_nim = poly_to_nim(bs[i]);
   ans[i] = nim_to_poly((Nimber(a_nim) / Nimber(b_nim)).val());
  }
  return ans;
 }
};
