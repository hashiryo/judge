// harness: 各 algos/*.hpp が定義する Solver::run(queries) を計測する。
// I/O 変換は計測外で行い、 純粋なアルゴリズム時間のみを ALGO_TIME_NS に含める。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);

 // 入力 parse (計測外)
 int T;
 std::cin >> T;
 std::vector<std::pair<u64, int>> queries(T);
 for (auto& [A, k] : queries) std::cin >> A >> k;

 std::vector<u64> result;
 uint64_t best_ns = ~uint64_t(0);
 for (int rep = 0; rep < 1; ++rep) {
  auto t0 = chrono::steady_clock::now();
  result = Solver::run(queries);
  auto t1 = chrono::steady_clock::now();
  auto ns = (uint64_t) chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if (ns < best_ns) best_ns = ns;
 }
 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long) best_ns);

 // 出力 (計測外)
 for (auto& v : result) std::cout << v << '\n';
 return 0;
}
