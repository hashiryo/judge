#pragma once
#include "_common.hpp"
// 2/3/5/7 wheel sieve (period = 210, coprime residue 48 個):
//   1, 11, 13, 17, ..., 209 mod 210 だけを bit で持つ。
//   メモリ N · 48/210 / 8 = N · 1/35 byte (wheel30 の N/30 byte の ~85%、
//   bit_sieve の N/16 の ~46%)。
//   篩内ループも 48/210 ≈ 23% の値だけ touch するので wheel30 より ~12% 削減。
#pragma GCC optimize("O3,unroll-loops")
struct Solver {
 static constexpr u32 PERIOD = 210;
 static constexpr u32 COPRIME = 48;
 static constexpr auto RES = []() {
  std::array<u8, COPRIME> r{};
  u32 idx = 0;
  for (u32 i = 1; i < PERIOD; i += 2) if ((i % 2 && i % 3 && i % 5 && i % 7)) r[idx++] = (u8) i;
  return r;
 }();
 static constexpr auto GAP = []() {
  std::array<u8, COPRIME> g{};
  for (u32 i = 0; i < COPRIME; ++i) {
   u32 next = (i + 1 < COPRIME) ? RES[i + 1] : (RES[0] + PERIOD);
   g[i] = (u8)(next - RES[i]);
  }
  return g;
 }();
 static constexpr auto STATE_FLOOR = []() {
  std::array<u8, PERIOD> s{};
  u8 idx = 0;
  for (u32 i = 0; i < PERIOD; ++i) {
   s[i] = idx;
   if ((i % 2 && i % 3 && i % 5 && i % 7)) ++idx;
  }
  return s;
 }();
 static constexpr u32 to_ord(u32 x) { return (x / PERIOD) * COPRIME + STATE_FLOOR[x % PERIOD]; }
 static constexpr u32 to_val(u32 x) { return (x / COPRIME) * PERIOD + RES[x % COPRIME]; }

 static std::pair<u32, std::vector<u32>> run(u32 N, u32 A, u32 B) {

  const u32 num_ords = to_ord(N) + COPRIME;
  const u32 num_words = (num_ords + 63) / 64;
  std::vector<u64> bits(num_words, ~u64(0));
  bits[0] &= ~u64(1);  // value 1 は素数でない

  // RES の中で値 ≥ 11 (= RES[1]) から開始。 sp = 1 で p = 11
  for (u32 sp = 1, p = RES[1] /* = 11 */;; ) {
   if (u64(p) * p > N) break;
   if (bits[to_ord(p) >> 6] & (u64(1) << (to_ord(p) & 63))) {
    u32 sk = sp;
    u64 m = u64(p) * p;
    while (m <= N) {
     u32 idx = to_ord(u32(m));
     bits[idx >> 6] &= ~(u64(1) << (idx & 63));
     u32 g = GAP[sk];
     if (++sk == COPRIME) sk = 0;
     m += u64(p) * g;
    }
   }
   u32 g = GAP[sp];
   if (++sp == COPRIME) sp = 0;
   p += g;
  }

  std::vector<u32> primes;
  primes.reserve(N / 10 + 16);
  if (N >= 2) primes.push_back(2);
  if (N >= 3) primes.push_back(3);
  if (N >= 5) primes.push_back(5);
  if (N >= 7) primes.push_back(7);
  for (u32 i = 1; i < num_ords; ++i) {
   if (bits[i >> 6] & (u64(1) << (i & 63))) {
    u32 v = to_val(i);
    if (v > N) break;
    primes.push_back(v);
   }
  }

  u32 cnt = u32(primes.size());
  u32 X = (cnt < B) ? 0 : (cnt - B + A - 1) / A;
  std::vector<u32> selected(X);
  for (u32 k = 0; k < X; ++k) selected[k] = primes[B + k * A];
  return {cnt, std::move(selected)};
 }
};
