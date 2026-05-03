#pragma once
#include "_common.hpp"
// =============================================================================
// yosupo Rational Approximation を cuixiang19 流の inline mediant ループで
// 実装したバリアント (xiaoziyao の SBT<T>::binary_search の代替実験)。
//
// 元: cuixiang19 の stern_brocot_tree fastest 提出 (https://judge.yosupo.jp/submission/327748)
// と同じく (lx, ly, rx, ry) を持って mediant 計算で進む。
// 各ステップで進める最大 k を直接除算で求めて一気に進める:
//   - x/y を跨がない最大 k:   k*(rx*y - ry*x) ≤ ly*x - lx*y
//   - 分母/分子 ≤ n を保つ:   lx + k*rx ≤ n, ly + k*ry ≤ n
// =============================================================================
struct Solver {
 static std::tuple<i64, i64, i64, i64> solve_one(i64 n, i64 x, i64 y) {
  i64 lx = 0, ly = 1, rx = 1, ry = 0;
  while (true) {
   // 既に左/右の片方が x/y と一致している場合 (前ステップで k_cross が exact 着地)
   if (lx > 0 && lx * y == ly * x) { rx = lx; ry = ly; break; }
   if (ry > 0 && rx * y == ry * x) { lx = rx; ly = ry; break; }
   i64 mx = lx + rx, my = ly + ry;
   if (mx > n || my > n) break;
   if (mx * y == my * x) { lx = rx = mx; ly = ry = my; break; }
   if (mx * y < my * x) {
    // 左を右に寄せる (R 方向): lx += k*rx, ly += k*ry
    // x/y を超えない最大 k: k*(rx*y - ry*x) ≤ ly*x - lx*y  (rx*y - ry*x > 0)
    i64 denom = rx * y - ry * x;
    i64 k_cross = denom > 0 ? (ly * x - lx * y) / denom : (i64) 4e18;
    i64 k_n = (n - lx) / rx;
    if (ry > 0) k_n = std::min(k_n, (n - ly) / ry);
    i64 k = std::min(k_cross, k_n);
    if (k <= 0) break;
    lx += k * rx;
    ly += k * ry;
   } else {
    // 右を左に寄せる (L 方向): rx += k*lx, ry += k*ly
    // x/y を下回らない最大 k: k*(ly*x - lx*y) ≤ rx*y - ry*x  (ly*x - lx*y > 0)
    i64 denom = ly * x - lx * y;
    i64 k_cross = denom > 0 ? (rx * y - ry * x) / denom : (i64) 4e18;
    i64 k_n = (n - ry) / ly;
    if (lx > 0) k_n = std::min(k_n, (n - rx) / lx);
    i64 k = std::min(k_cross, k_n);
    if (k <= 0) break;
    rx += k * lx;
    ry += k * ly;
   }
  }
  return {lx, ly, rx, ry};
 }

 static std::vector<std::tuple<i64, i64, i64, i64>> run(const std::vector<std::tuple<i64, i64, i64>>& queries) {
  std::vector<std::tuple<i64, i64, i64, i64>> ans;
  ans.reserve(queries.size());
  for (auto& [n, x, y] : queries) ans.push_back(solve_one(n, x, y));
  return ans;
 }
};
