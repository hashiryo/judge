// harness: 各 algos/*.hpp が定義する Solver::run(queries) を計測する。
// I/O 変換は計測外で行う。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/fastest.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);

 int T;
 std::cin >> T;
 std::vector<std::tuple<i64, i64, i64>> queries(T);
 for (auto& [n, x, y] : queries) std::cin >> n >> x >> y;

 std::vector<std::tuple<i64, i64, i64, i64>> result;
 uint64_t best_ns = ~uint64_t(0);
 for (int rep = 0; rep < 1; ++rep) {
  auto t0 = chrono::steady_clock::now();
  result = Solver::run(queries);
  auto t1 = chrono::steady_clock::now();
  auto ns = (uint64_t) chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if (ns < best_ns) best_ns = ns;
 }
 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long) best_ns);

 for (auto& [a, b, c, d] : result) std::cout << a << ' ' << b << ' ' << c << ' ' << d << '\n';
 return 0;
}
