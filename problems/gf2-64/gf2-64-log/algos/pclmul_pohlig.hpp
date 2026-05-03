#pragma once
// Pohlig-Hellman + BSGS で離散対数 log_2(x) in F_{2^64}^*。
//
// Group order N = 2^64 - 1 = 3 · 5 · 17 · 257 · 641 · 65537 · 6700417 (smooth)。
// 各素因数 p について:
//   1) x_p = x^{N/p}  (= subgroup of order p に射影)
//   2) g_p = g^{N/p}
//   3) BSGS で log_{g_p}(x_p) ∈ [0, p) を求める
//   4) CRT で全部合成 → log_g(x) mod N
//
// 既存 pclmul.hpp との違い:
//   - BSGS の baby-step テーブルは g に依存し x に依存しない → init で 1 回だけ構築
//   - flat hash map (open addressing, linear probing) で std::unordered_map より高速
//   - byte_window pow で x^{N/p} 計算を高速化
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_log_pohlig {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::reduce;

PCLMUL_FN u64 spread_bits(u32 a) {
#if defined(__BMI2__) && !defined(USE_SIMDE)
 return _pdep_u64(u64(a), 0x5555555555555555ull);
#else
 u64 x = a;
 x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x = (x | (x <<  8)) & 0x00FF00FF00FF00FFull;
 x = (x | (x <<  4)) & 0x0F0F0F0F0F0F0F0Full;
 x = (x | (x <<  2)) & 0x3333333333333333ull;
 x = (x | (x <<  1)) & 0x5555555555555555ull;
 return x;
#endif
}
PCLMUL_FN u64 sq(u64 a) {
 const u64 lo = spread_bits(u32(a));
 const u64 hi = spread_bits(u32(a >> 32));
 return reduce(__m128i{(long long) lo, (long long) hi});
}

inline u64 FROB4_BYTE[8][256];
inline bool inited = false;

PCLMUL_FN u64 frob4(u64 a) {
 return FROB4_BYTE[0][u8(a)]       ^ FROB4_BYTE[1][u8(a >>  8)]
      ^ FROB4_BYTE[2][u8(a >> 16)] ^ FROB4_BYTE[3][u8(a >> 24)]
      ^ FROB4_BYTE[4][u8(a >> 32)] ^ FROB4_BYTE[5][u8(a >> 40)]
      ^ FROB4_BYTE[6][u8(a >> 48)] ^ FROB4_BYTE[7][u8(a >> 56)];
}

// byte_window pow (= pclmul_byte_window.hpp と同じ)
PCLMUL_FN u64 pow_bw(u64 a, u64 e) {
 if (e == 0) return 1;
 u64 T[16];
 T[0] = 1; T[1] = a;
 #pragma GCC unroll 14
 for (int i = 2; i < 16; ++i) T[i] = mul(T[i-1], a);
 int top = 15;
 while (top > 0 && ((e >> (4 * top)) & 0xF) == 0) --top;
 u64 acc = T[(e >> (4 * top)) & 0xF];
 for (int i = top - 1; i >= 0; --i) {
  acc = frob4(acc);
  unsigned chunk = unsigned((e >> (4 * i)) & 0xF);
  if (chunk) acc = mul(acc, T[chunk]);
 }
 return acc;
}

// Open addressing hash map (linear probing) for BSGS
struct BSGSTable {
 std::vector<u64> keys;
 std::vector<u32> values;
 u64 mask;  // = capacity - 1, must be power of 2

 PCLMUL_RUN void build(u64 base, u64 m) {
  // capacity = next pow2 >= 4*m for low load factor
  u64 cap = 8;
  while (cap < 4 * m) cap *= 2;
  keys.assign(cap, ~u64(0));
  values.assign(cap, 0);
  mask = cap - 1;
  u64 cur = 1;
  for (u64 j = 0; j < m; ++j) {
   u64 h = (cur * 0x9E3779B97F4A7C15ull) & mask;
   while (keys[h] != ~u64(0)) h = (h + 1) & mask;
   keys[h] = cur;
   values[h] = u32(j);
   cur = mul(cur, base);
  }
 }

 PCLMUL_RUN u64 lookup(u64 key) const {
  u64 h = (key * 0x9E3779B97F4A7C15ull) & mask;
  while (keys[h] != ~u64(0)) {
   if (keys[h] == key) return values[h];
   h = (h + 1) & mask;
  }
  return ~u64(0);  // not found
 }
};

constexpr u64 N_ORDER = ~u64(0);  // 2^64 - 1
constexpr u64 G = 2;
constexpr u64 ORDER_PRIMES[7] = {3, 5, 17, 257, 641, 65537, 6700417};

struct BSGSCtx {
 u64 base;        // = g^{N/p}
 u64 inv_base_m;  // = base^{-m} (for giant step)
 u64 m;           // ceil(sqrt(p))
 u64 q;           // = p
 BSGSTable table;
};

inline BSGSCtx ctxs[7];

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 // frob4 byte table for pow_bw
 {
  u64 col[64];
  for (int j = 0; j < 64; ++j) {
   u64 v = u64(1) << j;
   for (int k = 0; k < 4; ++k) v = sq(v);
   col[j] = v;
  }
  for (int p = 0; p < 8; ++p) {
   for (int b = 0; b < 256; ++b) {
    u64 v = 0;
    for (int bit = 0; bit < 8; ++bit) {
     if ((b >> bit) & 1) v ^= col[p * 8 + bit];
    }
    FROB4_BYTE[p][b] = v;
   }
  }
 }
 // BSGS context for each prime
 for (int i = 0; i < 7; ++i) {
  const u64 p = ORDER_PRIMES[i];
  const u64 N_div_p = N_ORDER / p;
  ctxs[i].q = p;
  ctxs[i].base = pow_bw(G, N_div_p);
  u64 m = 1;
  while (m * m < p) ++m;
  ctxs[i].m = m;
  ctxs[i].inv_base_m = pow_bw(ctxs[i].base, p - m);  // = base^{-m} in subgroup of order p
  ctxs[i].table.build(ctxs[i].base, m);
 }
}

// Solve g_sub^k = target in subgroup of order q. k ∈ [0, q).
PCLMUL_FN u64 bsgs_solve(int prime_idx, u64 target) {
 const auto& ctx = ctxs[prime_idx];
 u64 t = target;
 for (u64 i = 0; i < ctx.m; ++i) {
  u64 j = ctx.table.lookup(t);
  if (j != ~u64(0)) {
   u64 res = i * ctx.m + j;
   if (res < ctx.q) return res;
  }
  t = mul(t, ctx.inv_base_m);
 }
 return ~u64(0);  // failure (shouldn't happen for valid target)
}

// Modular inverse using extended Euclidean (small numbers)
constexpr u64 modinv(u64 a, u64 m) {
 long long r0 = (long long) m, r1 = (long long) a;
 long long s0 = 0, s1 = 1;
 while (r1 != 0) {
  long long q = r0 / r1;
  long long r2 = r0 - q * r1; r0 = r1; r1 = r2;
  long long s2 = s0 - q * s1; s0 = s1; s1 = s2;
 }
 if (s0 < 0) s0 += (long long) m;
 return (u64) s0;
}

PCLMUL_FN u64 log_g(u64 x) {
 // Pohlig-Hellman: r ≡ r_i (mod p_i) で再構成
 // CRT: r = sum (r_i * M_i * y_i)  where M_i = ∏_{j≠i} p_j, y_i = M_i^{-1} mod p_i
 // 簡便な漸進 CRT: r mod (p_1 ... p_k) を更新
 u64 result = 0;
 u64 mod = 1;
 u64 r_arr[7];
 #pragma GCC unroll 7
 for (int i = 0; i < 7; ++i) {
  const u64 p = ORDER_PRIMES[i];
  const u64 x_sub = pow_bw(x, N_ORDER / p);
  r_arr[i] = bsgs_solve(i, x_sub);
 }
 // 漸進 CRT
 for (int i = 0; i < 7; ++i) {
  const u64 p = ORDER_PRIMES[i];
  u64 r_i = r_arr[i];
  u64 cur_mod_p = result % p;
  u64 diff = (r_i >= cur_mod_p) ? (r_i - cur_mod_p) : (p - (cur_mod_p - r_i));
  u64 mod_mod_p = mod % p;
  u64 inv = modinv(mod_mod_p, p);
  u64 t = (__uint128_t(diff) * inv) % p;
  result += mod * t;
  mod *= p;
 }
 return result;
}

} // namespace gf2_64_log_pohlig

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& xs) {
  using gf2_64_log_pohlig::log_g;
  using gf2_64_log_pohlig::init_tables;
  init_tables();
  vector<u64> ans(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) ans[i] = log_g(xs[i]);
  return ans;
 }
};
