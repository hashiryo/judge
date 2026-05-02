#pragma once
// 自前ライブラリ (judge/lib/mylib) の Nimber を呼ぶラッパー。
// 16-bit subfield 分解 + log/exp テーブル (65536 entry) で
// 4 つの 16-bit nim 積を組み合わせて 64-bit nim 積を実現。
#include "_common.hpp"
#include "mylib/algebra/Nimber.hpp"

struct NimProduct {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  Nimber::init(); // log/exp テーブル初期化 (1 度だけ)
  const size_t T = as.size();
  vector<u64> ans(T);
  for (size_t i = 0; i < T; ++i) {
   ans[i] = (Nimber(as[i]) * Nimber(bs[i])).val();
  }
  return ans;
 }
};
