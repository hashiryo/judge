#pragma once
// PCLMUL ベース版 + Pohlig-Hellman + BSGS。
// reference.hpp と同じ構造だが mul/pow に PCLMUL を使う。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#define PCLMUL_TARGET [[gnu::target("pclmul")]]
#else
#define PCLMUL_TARGET
#endif

namespace gf2_64_pcl_log {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::pow;

constexpr std::array<u64, 7> ORDER_PRIMES = {3, 5, 17, 257, 641, 65537, 6700417};

PCLMUL_TARGET inline u64 bsgs_subgroup(u64 base, u64 target, u64 q) {
 u64 m = 1;
 while (m * m < q) ++m;
 std::unordered_map<u64, u64> table;
 table.reserve(m * 2);
 u64 cur = 1;
 for (u64 j = 0; j < m; ++j) {
  table.try_emplace(cur, j);
  cur = mul(cur, base);
 }
 u64 inv_bm = pow(base, q - m);
 u64 t = target;
 for (u64 i = 0; i < m; ++i) {
  auto it = table.find(t);
  if (it != table.end()) {
   u64 res = i * m + it->second;
   if (res < q) return res;
  }
  t = mul(t, inv_bm);
 }
 return ~u64(0);
}

PCLMUL_TARGET inline u64 log_g(u64 x) {
 u64 N = GROUP_ORDER;
 u64 g = LOG_GENERATOR;
 u64 result = 0;
 u64 mod = 1;
 for (u64 q : ORDER_PRIMES) {
  u64 g_sub = pow(g, N / q);
  u64 x_sub = pow(x, N / q);
  u64 r = bsgs_subgroup(g_sub, x_sub, q);
  u64 mm = mod % q, e = q - 2;
  u64 inv_mod_q = 1;
  while (e) { if (e & 1) inv_mod_q = (__uint128_t) inv_mod_q * mm % q; mm = (__uint128_t) mm * mm % q; e >>= 1; }
  u64 diff = (r >= result % q) ? (r - result % q) : (q - (result % q - r));
  u64 t_ = (__uint128_t) diff * inv_mod_q % q;
  result += mod * t_;
  mod *= q;
 }
 return result;
}
}

struct GF2_64Op {
 PCLMUL_TARGET static vector<u64> run(const vector<u64>& xs) {
  vector<u64> ans(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) ans[i] = gf2_64_pcl_log::log_g(xs[i]);
  return ans;
 }
};
