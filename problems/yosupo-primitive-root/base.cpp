// harness: 各 algos/*.hpp が定義する struct PrimitiveRoot::run(qs) を計測する。
// yosupo "Primitive Root" 形式の I/O。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 int Q;
 cin >> Q;
 vector<u64> qs(Q);
 for(auto& x: qs) cin >> x;

 constexpr int REPEAT= 1;
 uint64_t best_ns= ~uint64_t(0);
 vector<u64> result;

 for(int rep= 0; rep < REPEAT; ++rep) {
  auto t0= chrono::steady_clock::now();
  auto r= PrimitiveRoot::run(qs);
  auto t1= chrono::steady_clock::now();
  result= std::move(r);
  auto ns= (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if(ns < best_ns) best_ns= ns;
 }

 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
 for(auto g: result) cout << g << '\n';
 return 0;
}
