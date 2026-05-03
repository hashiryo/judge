// harness: 各 algos/*.hpp が定義する struct Solver::run(input) を計測する。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 std::string input;
 {
  std::ostringstream oss;
  oss << std::cin.rdbuf();
  input = std::move(oss).str();
 }
 std::string output;
 uint64_t best_ns = ~uint64_t(0);
 for (int rep = 0; rep < 1; ++rep) {
  auto t0 = chrono::steady_clock::now();
  output = Solver::run(input);
  auto t1 = chrono::steady_clock::now();
  auto ns = (uint64_t) chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if (ns < best_ns) best_ns = ns;
 }
 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long) best_ns);
 std::cout << output;
 return 0;
}
