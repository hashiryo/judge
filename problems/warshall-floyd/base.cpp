// harness: 各 algos/*.hpp が定義する struct WF::run() を計測する。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 int V_;
 u64 W_, seed_;
 cin >> V_ >> W_ >> seed_;
 const int Vp= (V_ + 7) & ~7;  // V を 8 の倍数に切り上げ

 // 距離行列を Vp x Vp で確保し、INF で初期化、対角は 0。
 // 非対角 (i, j) は LCG で生成した重み (1..W)。パディング行/列は INF のまま。
 vector<i32> d((size_t)Vp * Vp, WF_INF);
 {
  constexpr u64 LCG_MA= 6364136223846793005ULL;
  constexpr u64 LCG_MB= 1442695040888963407ULL;
  u64 s= seed_;
  for(int i= 0; i < V_; ++i) {
   d[(size_t)i * Vp + i]= 0;
   for(int j= 0; j < V_; ++j) {
    if(i == j) continue;
    s= s * LCG_MA + LCG_MB;
    d[(size_t)i * Vp + j]= (i32)((s >> 16) % W_) + 1;
   }
  }
 }

 constexpr int REPEAT= 1;
 uint64_t best_ns= ~uint64_t(0);
 u64 result_out= 0;

 // d 初期状態を保持して REPEAT に備える (今は REPEAT=1 なので意味は薄いが将来対応)。
 vector<i32> d0= d;

 for(int rep= 0; rep < REPEAT; ++rep) {
  d= d0;
  WF wf(V_, Vp, d.data());
  auto t0= chrono::steady_clock::now();
  wf.run();
  auto t1= chrono::steady_clock::now();
  u64 acc= 0;
  for(int i= 0; i < V_; ++i)
   for(int j= 0; j < V_; ++j) acc^= (u64)(u32)d[(size_t)i * Vp + j];
  result_out= acc;
  auto ns= (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
  if(ns < best_ns) best_ns= ns;
 }

 fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
 cout << result_out << '\n';
 return 0;
}
