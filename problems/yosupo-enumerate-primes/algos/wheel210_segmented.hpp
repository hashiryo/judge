#pragma once
#include "_common.hpp"
// 2/3/5/7 wheel + segmented sieve:
//   block 単位で篩を打って cache 内で完結。 大 N で wheel210 から更に速化。
#pragma GCC optimize("O3,unroll-loops")
struct Solver {
 static constexpr u32 PERIOD = 210;
 static constexpr u32 COPRIME = 48;
 static constexpr auto RES = []() {
  std::array<u8, COPRIME> r{};
  u32 idx = 0;
  for (u32 i = 1; i < PERIOD; i += 2) if (i % 2 && i % 3 && i % 5 && i % 7) r[idx++] = (u8) i;
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
   if (i % 2 && i % 3 && i % 5 && i % 7) ++idx;
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

  // Step 1: small primes ≤ √N
  u32 sqrt_N = u32(std::sqrt(double(N))) + 2;
  while (u64(sqrt_N) * sqrt_N > N) --sqrt_N;
  std::vector<u32> sp;
  std::vector<u8> sp_state;
  {
   std::vector<u8> tmp(sqrt_N + 1, 1);
   tmp[0] = tmp[1] = 0;
   for (u32 i = 2; i * i <= sqrt_N; ++i) if (tmp[i]) {
    for (u32 j = i*i; j <= sqrt_N; j += i) tmp[j] = 0;
   }
   for (u32 i = 11; i <= sqrt_N; ++i) {
    if (tmp[i]) {
     sp.push_back(i);
     // find state
     u32 r = i % PERIOD;
     u8 s = 0;
     for (u8 j = 0; j < COPRIME; ++j) if (RES[j] == r) { s = j; break; }
     sp_state.push_back(s);
    }
   }
  }
  // 各 sp[i] の next position と state
  std::vector<u64> next_m(sp.size());
  std::vector<u8> next_s(sp.size());
  for (size_t i = 0; i < sp.size(); ++i) {
   next_m[i] = u64(sp[i]) * sp[i];
   next_s[i] = sp_state[i];  // state は k = p で開始
  }

  // Block: 2^22 values per block
  constexpr u32 BLOCK_VAL = 1u << 22;
  const u32 BLOCK_ORDS = to_ord(BLOCK_VAL) + COPRIME;
  const u32 BLOCK_WORDS = (BLOCK_ORDS + 63) / 64;

  std::vector<u32> primes;
  primes.reserve(N / 10 + 16);
  if (N >= 2) primes.push_back(2);
  if (N >= 3) primes.push_back(3);
  if (N >= 5) primes.push_back(5);
  if (N >= 7) primes.push_back(7);

  std::vector<u64> bits(BLOCK_WORDS);
  for (u64 lo = 0; lo <= N; lo += BLOCK_VAL) {
   u64 hi = std::min<u64>(lo + BLOCK_VAL, u64(N) + 1);
   std::fill(bits.begin(), bits.end(), ~u64(0));
   if (lo == 0) bits[0] &= ~u64(1);  // value 1 は素数でない

   for (size_t i = 0; i < sp.size(); ++i) {
    u64 m = next_m[i];
    u8 s = next_s[i];
    u32 p = sp[i];
    while (m < hi) {
     u32 idx = to_ord(u32(m)) - to_ord(u32(lo));
     bits[idx >> 6] &= ~(u64(1) << (idx & 63));
     u32 g = GAP[s];
     if (++s == COPRIME) s = 0;
     m += u64(p) * g;
    }
    next_m[i] = m;
    next_s[i] = s;
   }

   const u32 ord_hi = (hi == u64(N) + 1) ? to_ord(N) + 1 : to_ord(u32(hi));
   const u32 ord_lo = to_ord(u32(lo));
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
  out << cnt << ' ' << X << '\n';
  for (u32 k = 0; k < X; ++k) out << primes[B + k * A] << " \n"[k + 1 == X];
  if (X == 0) out << '\n';
  return std::move(out).str();
 }
};
