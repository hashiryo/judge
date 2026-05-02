#pragma once
// 塔分解アプローチ: poly → nim → Nimber::sqrt → poly。
// Nimber.hpp の sqrt は 16-bit subfield 内で sqrt + 上位 bit の補正を組み合わせた
// 専用の高速実装。pclmul.hpp は a^(2^63) 経由なので tower のほうが速いはず。
#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/basis_change.hpp"
#include "mylib/algebra/Nimber.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_basis::poly_to_nim;
  using gf2_64_basis::nim_to_poly;
  Nimber::init();
  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   u64 a_nim = poly_to_nim(as[i]);
   ans[i] = nim_to_poly(Nimber(a_nim).sqrt().val());
  }
  return ans;
 }
};
