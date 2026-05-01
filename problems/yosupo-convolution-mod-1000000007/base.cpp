// harness: 各 algos/*.hpp が定義する struct Conv::run(const vector<u32>&, const vector<u32>&)
// を計測する。yosupo の "Convolution mod 10^9+7" 形式の I/O。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 int N, M;
 cin >> N >> M;
 vector<u32> a(N), b(M);
 for(auto& x: a) cin >> x;
 for(auto& x: b) cin >> x;

 constexpr int REPEAT= 1;
 uint64_t best_ns= ~uint64_t(0);
 vector<u32> result;

 for(int rep= 0; rep < REPEAT; ++rep) {
  auto t0= chrono::steady_clock::now();
  auto r= Conv::run(a, b);
  auto t1= chrono::steady_clock::now();
  result= std::move(r);
  auto ns= (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if(ns < best_ns) best_ns= ns;
 }

 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
 for(size_t i= 0; i < result.size(); ++i) {
  if(i) cout << ' ';
  cout << result[i];
 }
 cout << '\n';
 return 0;
}
