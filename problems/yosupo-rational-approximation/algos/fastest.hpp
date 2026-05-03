#pragma once
#include "_common.hpp"
// =============================================================================
// yosupo Rational Approximation の最速提出 (xiaoziyao, https://judge.yosupo.jp/submission/290530)
// を std::string ベースの harness に移植したもの。
//
// アルゴリズム:
//   Stern-Brocot 木上の二分探索 (SBT::binary_search):
//     根 (1/1) から始め、各ステップで (a/b, c/d) を保持。
//     mediant (a+c)/(b+d) が条件を満たす場合と満たさない場合で再帰的に分岐。
//     denominator ≤ n の制約下で最良近似を探索。
//
//   各クエリ: O(log min(n, max(x,y))) 程度。
// =============================================================================
struct Solver {
 template<class T>
 struct SBT {
  // n 以下の denominator/numerator で条件 f を満たす最良 fraction の範囲を返す
  // f(a, b) は単調 (true → false の境界を探す)
  template<class F>
  static std::tuple<T, T, T, T> binary_search(T n, F f) {
   auto dfs = [&](auto& self, bool side, T& a, T& b, T c, T d) -> bool {
    auto ok = [&]() { return a + c <= n && b + d <= n && f(a + c, b + d) == side; };
    if (!ok()) return false;
    if (!self(self, side, a, b, c + c, d + d) || ok()) { a += c; b += d; }
    return true;
   };
   T a = 0, b = 1, c = 1, d = 0;
   while (dfs(dfs, true, a, b, c, d) || dfs(dfs, false, c, d, a, b)) {}
   return {a, b, c, d};
  }
 };

 static std::string run(const std::string& input) {
  std::istringstream in(input);
  std::ostringstream out;
  int t; in >> t;
  while (t--) {
   i64 n, x, y;
   in >> n >> x >> y;
   auto [a, b, c, d] = SBT<i64>::binary_search(n, [&](i64 a_, i64 b_) { return a_ * y <= b_ * x; });
   if (a * y == b * x) { c = a; d = b; }
   out << a << ' ' << b << ' ' << c << ' ' << d << '\n';
  }
  return std::move(out).str();
 }
};
