#pragma once
#include "_common.hpp"
// Bit-packed sieve + 「奇数のみ」表現:
//   odd_idx(2i+1) = i, ビット 1 = prime。 メモリは N/16 byte で naive の 1/16。
//   cache 効率が良くなり naive より大幅高速化。
struct Solver {
 static std::string run(const std::string& input) {
  std::istringstream in(input);
  std::ostringstream out;
  u32 N, A, B;
  in >> N >> A >> B;

  // bit i represents the number 2i+1 (odd values only)
  // We need indexes for i in [0, (N+1)/2], so words = (M + 63) / 64 where M = (N+1)/2 + 1
  const u32 M = (N - 1) / 2 + 1;  // odd numbers 3, 5, ..., are at indices 1, 2, ...
  // Actually let's use: index i represents value 2i+1. So i=0→1, i=1→3, i=2→5...
  // Want index for v ≤ N: max i = N/2 if N is odd, else (N-1)/2. Use (N+1)/2 to be safe.
  const u32 max_i = N / 2 + 1;
  const u32 words = (max_i + 64) / 64;
  std::vector<u64> bits(words, ~u64(0));
  // bit 0 (= value 1) is not prime
  bits[0] &= ~u64(1);

  // Sieve: for each odd p starting from 3, p^2 step 2p
  for (u32 p = 3; u64(p) * p <= N; p += 2) {
   u32 pi = p / 2;
   if (bits[pi >> 6] & (u64(1) << (pi & 63))) {
    // mark multiples of p starting from p*p, step 2p (only odd multiples)
    for (u64 j = u64(p) * p; j <= N; j += 2 * p) {
     u32 ji = u32(j / 2);
     bits[ji >> 6] &= ~(u64(1) << (ji & 63));
    }
   }
  }

  // 全 prime を抽出。 2 だけ別扱い。
  std::vector<u32> primes;
  primes.reserve(N / 10 + 16);
  if (N >= 2) primes.push_back(2);
  for (u32 i = 1; 2 * i + 1 <= N; ++i) {
   if (bits[i >> 6] & (u64(1) << (i & 63))) primes.push_back(2 * i + 1);
  }

  u32 cnt = u32(primes.size());
  u32 X = (cnt < B) ? 0 : (cnt - B + A - 1) / A;
  out << cnt << ' ' << X << '\n';
  for (u32 k = 0; k < X; ++k) out << primes[B + k * A] << " \n"[k + 1 == X];
  if (X == 0) out << '\n';
  return std::move(out).str();
 }
};
