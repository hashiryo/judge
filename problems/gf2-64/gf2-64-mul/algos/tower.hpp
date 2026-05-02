#pragma once
// 塔分解アプローチ:
//   poly 基底 u64 → 基底変換で nim 自然表現 → Nimber::mul (mylib, F_{2^16}^4 塔)
//   → 基底変換で poly 基底 u64 に戻す
//
// nim 自然表現側で計算するため、以下の利点がある:
//   - mylib の Nimber は 16-bit subfield + log/exp 256 KB テーブルで高速
//   - PCLMUL を使わないので PCLMUL 非対応環境でも動く
// 一方、poly canonical なので前後に **基底変換 8 byte lookup × 2 direction = 16 lookup**
// のオーバーヘッドが乗る。
//
// pclmul.hpp との比較ポイント:
//   - x86 (PCLMUL あり): pclmul.hpp が圧勝するはず
//   - ARM/PCLMUL なし: tower.hpp が pclmul.hpp (SIMDe 経由) に勝つ可能性

#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/basis_change.hpp"
#include "mylib/algebra/Nimber.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_basis::poly_to_nim;
  using gf2_64_basis::nim_to_poly;
  Nimber::init();  // 16-bit subfield log/exp テーブル構築 (1 度だけ)
  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   u64 a_nim = poly_to_nim(as[i]);
   u64 b_nim = poly_to_nim(bs[i]);
   ans[i] = nim_to_poly((Nimber(a_nim) * Nimber(b_nim)).val());
  }
  return ans;
 }
};
