#pragma once
#include "_common.hpp"
// wheel30 + Segmented sieve:
//   篩配列を block (~256KB = L2 fit) ごとに区切って処理する。
//   各 block は (bit array of one block) で完結、 cache 内で sieve 完了。
//
// 期待: 大 N で wheel30 から ~1.5-2× 高速化 (wheel30 は最終 N/30 byte が L3 spill)
#pragma GCC optimize("O3,unroll-loops")
struct Solver {
 static constexpr u32 PERIOD = 30;
 static constexpr u32 COPRIME = 8;
 static constexpr std::array<u8, COPRIME> RES = {1, 7, 11, 13, 17, 19, 23, 29};
 static constexpr std::array<u8, COPRIME> GAP = {6, 4, 2, 4, 2, 4, 6, 2};
 static constexpr auto STATE_FLOOR = []() {
  std::array<u8, PERIOD> s{};
  u8 idx = 0;
  for (u32 i = 0; i < PERIOD; ++i) {
   s[i] = idx;
   if (i % 2 && i % 3 && i % 5) ++idx;
  }
  return s;
 }();
 static constexpr u32 to_ord(u32 x) { return (x / PERIOD) * COPRIME + STATE_FLOOR[x % PERIOD]; }
 static constexpr u32 to_val(u32 x) { return (x / COPRIME) * PERIOD + RES[x % COPRIME]; }

 static std::pair<u32, std::vector<u32>> run(u32 N, u32 A, u32 B) {

  // Step 1: precompute small primes ≤ √N via simple bit sieve
  u32 sqrt_N = u32(std::sqrt(double(N))) + 2;
  while (u64(sqrt_N) * sqrt_N > N) --sqrt_N;
  std::vector<u32> sp;
  {
   std::vector<u8> tmp(sqrt_N + 1, 1);
   tmp[0] = tmp[1] = 0;
   for (u32 i = 2; i * i <= sqrt_N; ++i) if (tmp[i]) {
    for (u32 j = i*i; j <= sqrt_N; j += i) tmp[j] = 0;
   }
   for (u32 i = 7; i <= sqrt_N; ++i) {
    if (tmp[i]) sp.push_back(i);
   }
  }
  // 各 sp[i] の next-position と state を維持
  // next ≥ p^2 で、 coprime to 30 なもの
  std::vector<u64> next_m(sp.size());
  std::vector<u8> next_s(sp.size());
  for (size_t i = 0; i < sp.size(); ++i) {
   u64 m = u64(sp[i]) * sp[i];
   // state は k = p の residue index (m = p^2 のときも k = p なので)
   u32 r = sp[i] % PERIOD;
   u8 s = 0;
   for (u8 j = 0; j < COPRIME; ++j) if (RES[j] == r) { s = j; break; }
   next_m[i] = m;
   next_s[i] = s;
  }

  // Block size: bit array per block fits in L2 (~256KB). 256KB = 2M bits.
  // value range per block = 2M / (8/30) = 7.5M
  constexpr u32 BLOCK_VAL = 1u << 22;  // ~4M values per block
  const u32 BLOCK_ORDS = to_ord(BLOCK_VAL) + COPRIME;
  const u32 BLOCK_WORDS = (BLOCK_ORDS + 63) / 64;

  std::vector<u32> primes;
  primes.reserve(N / 10 + 16);
  if (N >= 2) primes.push_back(2);
  if (N >= 3) primes.push_back(3);
  if (N >= 5) primes.push_back(5);

  std::vector<u64> bits(BLOCK_WORDS);
  for (u32 lo = 0; lo <= N; lo += BLOCK_VAL) {
   u32 hi = std::min<u64>(u64(lo) + BLOCK_VAL, u64(N) + 1);
   // initialize block all-1
   std::fill(bits.begin(), bits.end(), ~u64(0));
   if (lo == 0) bits[0] &= ~u64(1);  // value 1 is not prime

   // Sieve: for each small prime, mark its multiples in [lo, hi)
   for (size_t i = 0; i < sp.size(); ++i) {
    u64 m = next_m[i];
    u8 s = next_s[i];
    u32 p = sp[i];
    while (m < hi) {
     // bit position relative to block: ord(m) - ord(lo)
     u32 idx = to_ord(u32(m)) - to_ord(lo);
     bits[idx >> 6] &= ~(u64(1) << (idx & 63));
     u32 g = GAP[s];
     s = (s + 1) & 7;
     m += u64(p) * g;
    }
    next_m[i] = m;
    next_s[i] = s;
   }

   // Collect primes from this block
   const u32 ord_hi = (hi == u64(N) + 1) ? to_ord(N) + 1 : to_ord(hi);
   const u32 ord_lo = to_ord(lo);
   for (u32 i = (lo == 0 ? 1 : 0); i < ord_hi - ord_lo; ++i) {
    if (bits[i >> 6] & (u64(1) << (i & 63))) {
     u32 v = to_val(ord_lo + i);
     if (v > N) break;
     primes.push_back(v);
    }
   }
  }

  u32 cnt = u32(primes.size());
  u32 X = (cnt < B) ? 0 : (cnt - B + A - 1) / A;
  std::vector<u32> selected(X);
  for (u32 k = 0; k < X; ++k) selected[k] = primes[B + k * A];
  return {cnt, std::move(selected)};
 }
};
