// harness: T 個の (a, b) を読み、a ⊗ b ∈ GF(2^64) を出力する。
//
// I/O:
//   T
//   a_0 b_0
//   ...
//   a_{T-1} b_{T-1}
// 出力: 各 a_i ⊗ b_i (= GF(2)[x] / (x^64+x^4+x^3+x+1) の積)
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/reference.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 int T;
 cin >> T;
 vector<u64> as(T), bs(T);
 for (int i = 0; i < T; ++i) cin >> as[i] >> bs[i];

 uint64_t best_ns = ~uint64_t(0);
 vector<u64> result;
 for (int rep = 0; rep < 1; ++rep) {
  auto t0 = chrono::steady_clock::now();
  auto r = GF2_64Op::run(as, bs);
  auto t1 = chrono::steady_clock::now();
  result = std::move(r);
  auto ns = (uint64_t) chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if (ns < best_ns) best_ns = ns;
 }
 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long) best_ns);
 for (auto x : result) cout << x << '\n';
 return 0;
}
