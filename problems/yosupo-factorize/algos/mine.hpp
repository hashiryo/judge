#pragma once
// 自前ライブラリ (judge/lib/mylib) の Factors を呼ぶラッパー。
// Pollard-Rho + Miller-Rabin 内蔵。

#include "_common.hpp"
// libc++ (macOS) には __lg が無いので fallback で提供 (mylib の Factors が依存)。
#if !defined(__GLIBCXX__)
namespace { template <class T> constexpr int __lg(T x) { return std::bit_width<std::make_unsigned_t<T>>(x) - 1; } }
#endif
#include "mylib/number_theory/Factors.hpp"

struct Factorize {
 static vector<vector<u64>> run(const vector<u64>& qs) {
  vector<vector<u64>> ans;
  ans.reserve(qs.size());
  for (auto x : qs) {
   vector<u64> fs;
   if (x > 1) {
    Factors f(x);
    for (auto [p, c] : f) {
     for (uint16_t k = 0; k < c; ++k) fs.push_back(p);
    }
    std::sort(fs.begin(), fs.end());
   }
   ans.push_back(std::move(fs));
  }
  return ans;
 }
};
