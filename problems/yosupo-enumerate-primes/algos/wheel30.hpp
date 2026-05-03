#pragma once
#include "_common.hpp"
// 2/3/5 wheel sieve (period = 30, coprime residue 8 個):
//   1, 7, 11, 13, 17, 19, 23, 29 mod 30 だけを bit で持つ。
//   メモリ N · 8/30 / 8 = N/30 byte (bit_sieve の N/16 の約半分)。
//
// 重要: inner loop の next-coprime 計算は precompute gap table で O(1)。
//   gap[s] = res[(s+1)%8] - res[s] mod 30 (with wrap)
//
// 期待: bit_sieve から ~1.3-1.5× 高速化 (メモリ削減 + 篩内ループ削減)
#pragma GCC optimize("O3,unroll-loops")
struct Solver {
 static constexpr u32 PERIOD = 30;
 static constexpr u32 COPRIME = 8;
 static constexpr std::array<u8, COPRIME> RES = {1, 7, 11, 13, 17, 19, 23, 29};
 // gap[s] = next residue - current residue (cyclic)
 // 1→7=6, 7→11=4, 11→13=2, 13→17=4, 17→19=2, 19→23=4, 23→29=6, 29→31=2
 static constexpr std::array<u8, COPRIME> GAP = {6, 4, 2, 4, 2, 4, 6, 2};
 // state[r mod 30] = (state_idx of residue ≥ r), or undefined if r is not coprime
 static constexpr auto STATE_OF_R = []() {
  std::array<u8, PERIOD> s{};
  for (u8 i = 0; i < COPRIME; ++i) s[RES[i]] = i;
  return s;
 }();
 // map value to ordinal (compressed bit position)
 static constexpr auto STATE_FLOOR = []() {
  // state_floor[r] = (state index of next coprime residue ≥ r)
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

 static std::string run(const std::string& input) {
  std::istringstream in(input);
  std::ostringstream out;
  u32 N, A, B;
  in >> N >> A >> B;

  const u32 num_ords = to_ord(N) + COPRIME;
  const u32 num_words = (num_ords + 63) / 64;
  std::vector<u64> bits(num_words, ~u64(0));
  bits[0] &= ~u64(1);  // value 1 は素数でない

  // 篩本体: p ∈ {7, 11, 13, ..., √N} を coprime gap table で巡回
  // 各 p について、 multiples p*k (k ≥ p, k coprime to 30) を消す
  for (u32 sp = 1, p = RES[1] /* = 7 */;; ) {
   if (u64(p) * p > N) break;
   if (bits[to_ord(p) >> 6] & (u64(1) << (to_ord(p) & 63))) {
    // p is prime, sieve multiples p*k for k coprime, starting k = p
    u32 sk = sp;  // k = p, residue index = sp
    u64 m = u64(p) * p;
    while (m <= N) {
     u32 idx = to_ord(u32(m));
     bits[idx >> 6] &= ~(u64(1) << (idx & 63));
     // advance k by GAP[sk]
     u32 g = GAP[sk];
     sk = (sk + 1) & 7;
     m += u64(p) * g;
    }
   }
   // advance p
   u32 g = GAP[sp];
   sp = (sp + 1) & 7;
   p += g;
  }

  std::vector<u32> primes;
  primes.reserve(N / 10 + 16);
  if (N >= 2) primes.push_back(2);
  if (N >= 3) primes.push_back(3);
  if (N >= 5) primes.push_back(5);
  for (u32 i = 1; i < num_ords; ++i) {
   if (bits[i >> 6] & (u64(1) << (i & 63))) {
    u32 v = to_val(i);
    if (v > N) break;
    primes.push_back(v);
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
