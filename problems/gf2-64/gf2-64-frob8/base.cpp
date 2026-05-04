// harness: T 個の u64 a_i を読み、 a_i^{2^8} = a_i^256 ∈ GF(2^64) を出力。
// I/O 変換は計測外、 GF2_64Op::run(as) のみ algo time に含める。

#ifndef ALGO_HPP
#define ALGO_HPP "algos/reference.hpp"
#endif
#include ALGO_HPP

signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 int T;
 cin >> T;
 vector<u64> as(T);
 for(auto& x: as) cin >> x;

 uint64_t best_ns= ~uint64_t(0);
 vector<u64> result;
 for(int rep= 0; rep < 1; ++rep) {
  auto t0= chrono::steady_clock::now();
  auto r= GF2_64Op::run(as);
  auto t1= chrono::steady_clock::now();
  result= std::move(r);
  auto ns= (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if(ns < best_ns) best_ns= ns;
 }
 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
 for(auto y: result) cout << y << '\n';
 return 0;
}
