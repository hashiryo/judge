// harness: 各 algos/*.hpp が定義する struct DLog::run(p, g, queries) を計測する。
// yosupo "Discrete Logarithm Fixed Mod" 形式の I/O。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 u32 p, g;
 int n;
 cin >> p >> g >> n;
 vector<u32> qs(n);
 for (auto& a : qs) cin >> a;

 uint64_t best_ns = ~uint64_t(0);
 vector<u32> result;
 for (int rep = 0; rep < 1; ++rep) {
  auto t0 = chrono::steady_clock::now();
  auto r = DLog::run(p, g, qs);
  auto t1 = chrono::steady_clock::now();
  result = std::move(r);
  auto ns = (uint64_t) chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if (ns < best_ns) best_ns = ns;
 }
 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long) best_ns);
 // u32(-1) を sentinel として "log なし" を表すため、int32_t 経由でキャスト
 for (auto x : result) cout << (int32_t) x << '\n';
 return 0;
}
