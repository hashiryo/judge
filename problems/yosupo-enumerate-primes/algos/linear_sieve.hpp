#pragma once
#include "_common.hpp"
// Linear sieve: 各合成数を smallest prime factor で 1 回だけ落とす O(N) 篩。
//   既知 prime list を順に走査し、 i * primes[j] の篩を打つ。
//   primes[j] | i になった時点で break することで重複を排除。
struct Solver {
 static std::string run(const std::string& input) {
  std::istringstream in(input);
  std::ostringstream out;
  u32 N, A, B;
  in >> N >> A >> B;

  std::vector<u32> spf(N + 1, 0);
  std::vector<u32> primes;
  primes.reserve(N / 10 + 16);
  for (u32 i = 2; i <= N; ++i) {
   if (spf[i] == 0) { spf[i] = i; primes.push_back(i); }
   for (u32 p : primes) {
    if (p > spf[i] || u64(p) * i > N) break;
    spf[p * i] = p;
   }
  }

  u32 cnt = u32(primes.size());
  u32 X = (cnt < B) ? 0 : (cnt - B + A - 1) / A;
  out << cnt << ' ' << X << '\n';
  for (u32 k = 0; k < X; ++k) out << primes[B + k * A] << " \n"[k + 1 == X];
  if (X == 0) out << '\n';
  return std::move(out).str();
 }
};
