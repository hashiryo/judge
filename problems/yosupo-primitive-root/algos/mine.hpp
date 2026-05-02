#pragma once
// 自前ライブラリ (judge/lib/mylib) の OrderFp::primitive_root を呼ぶラッパー。
// p-1 の Pollard-Rho 素因数分解 + witness 探索 (Mont 累乗) を内蔵。
#include "_common.hpp"
// libc++ (macOS) には __lg が無いので fallback で提供 (mylib が依存)。
#if !defined(__GLIBCXX__)
namespace { template <class T> constexpr int __lg(T x) { return std::bit_width<std::make_unsigned_t<T>>(x) - 1; } }
#endif
#include "mylib/number_theory/OrderFp.hpp"

struct PrimitiveRoot {
 static vector<u64> run(const vector<u64>& qs) {
  vector<u64> ans;
  ans.reserve(qs.size());
  for (auto p : qs) ans.push_back(OrderFp(p).primitive_root());
  return ans;
 }
};
