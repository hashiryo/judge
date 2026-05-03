// harness: 各 algos/*.hpp が定義する Solver::run(n, m, a, b) を計測する。
// 入出力フォーマット変換 (string ↔ vector<u64>) は計測外で行い、
// 純粋なアルゴリズム時間のみを ALGO_TIME_NS に含める。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);

 // 入力 parse (計測外)
 int n, m;
 std::cin >> n >> m;
 std::vector<u64> a(n), b(m);
 for (auto& x : a) std::cin >> x;
 for (auto& x : b) std::cin >> x;

 std::vector<u64> result;
 uint64_t best_ns = ~uint64_t(0);
 for (int rep = 0; rep < 1; ++rep) {
  auto t0 = chrono::steady_clock::now();
  result = Solver::run(n, m, a, b);
  auto t1 = chrono::steady_clock::now();
  auto ns = (uint64_t) chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if (ns < best_ns) best_ns = ns;
 }
 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long) best_ns);

 // 出力 (計測外)
 for (size_t k = 0; k < result.size(); ++k) {
  std::cout << result[k];
  std::cout << (k + 1 == result.size() ? '\n' : ' ');
 }
 return 0;
}
