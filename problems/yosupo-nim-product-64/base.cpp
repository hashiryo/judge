// harness: 各 algos/*.hpp が定義する struct NimProduct::run(as, bs) を計測する。
// yosupo "Nim Product (F_{2^64})" 形式の I/O。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 int T;
 cin >> T;
 vector<u64> as(T), bs(T);
 for (int i = 0; i < T; ++i) cin >> as[i] >> bs[i];

 constexpr int REPEAT = 1;
 uint64_t best_ns = ~uint64_t(0);
 vector<u64> result;

 for (int rep = 0; rep < REPEAT; ++rep) {
  auto t0 = chrono::steady_clock::now();
  auto r = NimProduct::run(as, bs);
  auto t1 = chrono::steady_clock::now();
  result = std::move(r);
  auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if (ns < best_ns) best_ns = ns;
 }

 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
 // 高速 stdout: 改行区切りで T 行
 for (auto x : result) cout << x << '\n';
 return 0;
}
