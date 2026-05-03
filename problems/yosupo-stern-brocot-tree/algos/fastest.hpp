#pragma once
#include "_common.hpp"
// =============================================================================
// yosupo Stern-Brocot Tree の最速提出 (cuixiang19, https://judge.yosupo.jp/submission/327748)
// を std::string ベースの I/O harness に移植したもの。
//
// アルゴリズム:
//   各クエリで Stern-Brocot 木上を mediant 計算で root から目標 rational に
//   到達する経路を辿る。LCA, ANCESTOR, RANGE 等は同じ traversal の派生。
//   全部が i64 整数の比較・加減算・除算で実装されているため極めて高速。
//
// 注: 元提出は luogu P1797 用 (= 同問題の中国版)、yosupo の入出力フォーマットと
//     同じ。
// =============================================================================
struct Solver {
 // 文字列パーサ: stringstream で operator>> を使う
 static std::string run(const std::string& input) {
  std::istringstream in(input);
  std::ostringstream out;
  int T;
  in >> T;
  while (T--) {
   solve_one(in, out);
  }
  return std::move(out).str();
 }

 static void solve_one(std::istream& in, std::ostream& out) {
  std::string op;
  in >> op;
  if (op == "DECODE_PATH") {
   int k; in >> k;
   i64 lx = 0, ly = 1, rx = 1, ry = 0;
   while (k--) {
    std::string dir; i64 x;
    in >> dir >> x;
    if (dir == "L") {
     rx += lx * x;
     ry += ly * x;
    } else {
     lx += rx * x;
     ly += ry * x;
    }
   }
   out << (lx + rx) << ' ' << (ly + ry) << '\n';
  } else if (op == "ENCODE_PATH") {
   i64 x, y; in >> x >> y;
   i64 lx = 0, ly = 1, rx = 1, ry = 0;
   std::vector<std::pair<char, i64>> ans;
   while (true) {
    i64 mx = lx + rx, my = ly + ry;
    if (mx == x && my == y) break;
    if (x * my < y * mx) {
     i64 k = (y * rx - x * ry - 1) / (x * ly - y * lx);
     ans.emplace_back('L', k);
     rx += lx * k;
     ry += ly * k;
    } else {
     i64 k = (x * ly - y * lx - 1) / (y * rx - x * ry);
     ans.emplace_back('R', k);
     lx += rx * k;
     ly += ry * k;
    }
   }
   out << (int) ans.size();
   for (auto [c, n] : ans) out << ' ' << c << ' ' << n;
   out << '\n';
  } else if (op == "LCA") {
   i64 ax, ay, bx, by;
   in >> ax >> ay >> bx >> by;
   if (ax * by > bx * ay) { std::swap(ax, bx); std::swap(ay, by); }
   i64 lx = 0, rx = 1, ly = 1, ry = 0;
   while (true) {
    i64 mx = lx + rx, my = ly + ry;
    if (ax * my <= ay * mx && mx * by <= my * bx) {
     out << mx << ' ' << my << '\n';
     break;
    }
    if (bx * my < by * mx) {
     i64 k = (by * rx - bx * ry - 1) / (bx * ly - by * lx);
     rx += lx * k;
     ry += ly * k;
    } else {
     i64 k = (ax * ly - ay * lx - 1) / (ay * rx - ax * ry);
     lx += rx * k;
     ly += ry * k;
    }
   }
  } else if (op == "ANCESTOR") {
   i64 z, x, y;
   in >> z >> x >> y;
   i64 lx = 0, ly = 1, rx = 1, ry = 0;
   std::vector<std::pair<char, i64>> ans;
   while (true) {
    i64 mx = lx + rx, my = ly + ry;
    if (mx == x && my == y) break;
    if (x * my < y * mx) {
     i64 k = (y * rx - x * ry - 1) / (x * ly - y * lx);
     ans.emplace_back('L', k);
     rx += lx * k;
     ry += ly * k;
    } else {
     i64 k = (x * ly - y * lx - 1) / (y * rx - x * ry);
     ans.emplace_back('R', k);
     lx += rx * k;
     ly += ry * k;
    }
   }
   lx = 0; ly = 1; rx = 1; ry = 0;
   for (auto [c, n] : ans) {
    if (c == 'L') {
     if (z <= n) { rx += lx * z; ry += ly * z; z = 0; break; }
     z -= n; rx += lx * n; ry += ly * n;
    } else {
     if (z <= n) { lx += rx * z; ly += ry * z; z = 0; break; }
     z -= n; lx += rx * n; ly += ry * n;
    }
   }
   if (z) { out << "-1\n"; return; }
   out << (lx + rx) << ' ' << (ly + ry) << '\n';
  } else {  // RANGE
   i64 x, y;
   in >> x >> y;
   i64 lx = 0, rx = 1, ly = 1, ry = 0;
   while (true) {
    i64 mx = lx + rx, my = ly + ry;
    if (mx == x && my == y) {
     out << lx << ' ' << ly << ' ' << rx << ' ' << ry << '\n';
     return;
    }
    if (x * my < y * mx) {
     i64 k = (y * rx - x * ry - 1) / (x * ly - y * lx);
     rx += lx * k;
     ry += ly * k;
    } else {
     i64 k = (x * ly - y * lx - 1) / (y * rx - x * ry);
     lx += rx * k;
     ly += ry * k;
    }
   }
  }
 }
};
