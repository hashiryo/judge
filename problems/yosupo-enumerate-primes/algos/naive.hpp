#pragma once
#include "_common.hpp"
// 素朴な Sieve of Eratosthenes (bool 配列)
struct Solver {
 static std::string run(const std::string& input) {
  std::istringstream in(input);
  std::ostringstream out;
  u32 N, A, B;
  in >> N >> A >> B;

  // is_composite[i] = true if i is composite
  std::vector<u8> is_comp(N + 1, 0);
  is_comp[0] = is_comp[1] = 1;
  for (u32 i = 2; u64(i) * i <= N; ++i) {
   if (!is_comp[i]) {
    for (u32 j = i * i; j <= N; j += i) is_comp[j] = 1;
   }
  }
  std::vector<u32> primes;
  primes.reserve(N / 10 + 16);
  for (u32 i = 2; i <= N; ++i) if (!is_comp[i]) primes.push_back(i);

  u32 cnt = u32(primes.size());
  u32 X = (cnt < B) ? 0 : (cnt - B + A - 1) / A;
  out << cnt << ' ' << X << '\n';
  for (u32 k = 0; k < X; ++k) out << primes[B + k * A] << " \n"[k + 1 == X];
  if (X == 0) out << '\n';
  return std::move(out).str();
 }
};
