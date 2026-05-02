#pragma once
// =============================================================================
// Source: yosupo "Primitive Root" 提出 336163.
//   https://judge.yosupo.jp/submission/336163
//   方式:
//     - Montgomery64 (mod < 2^62 前提) で剰余演算
//     - 5-base Miller-Rabin (ランダム base) で素数判定
//     - Pollard-Rho (Brent 風 batched gcd) で p-1 を素因数分解
//     - witness 探索: ランダム g を生成し、p-1 の各素因数 q について
//       g^((p-1)/q) != 1 を全て満たすまで試行
// 抽出方針:
//   - 元提出の I/O (cin/cout) は base.cpp に移譲
//   - Montgomery64 / is_prime / factorize / primitive_root をそのままコピー
//   - rd / rng は元実装と同じ random_device + mt19937_64
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include <random>
#include <map>
#include <numeric>

namespace yosupo_336163 {
using std::uint64_t;
using u128_local = __uint128_t;

struct Montgomery64 {
 uint64_t mod, mod_inv, one;
 Montgomery64(uint64_t mod): mod(mod), mod_inv(1) {
  for (int i = 0; i < 6; ++i) mod_inv *= 2 - mod * mod_inv;
  one = (u128_local(1) << 64) % mod;
 }
 uint64_t norm(uint64_t x) const { return x < x - mod ? x : x - mod; }
 uint64_t transform(uint64_t x) const { return ((u128_local) x << 64) % mod; }
 uint64_t reduce(u128_local x) const {
  uint64_t m = (uint64_t(x) * mod_inv * u128_local(mod)) >> 64;
  return norm((x >> 64) + mod - m);
 }
 uint64_t add(uint64_t x, uint64_t y) const { return x + y >= mod ? x + y - mod : x + y; }
 uint64_t sub(uint64_t x, uint64_t y) const { return x >= y ? x - y : x + mod - y; }
 uint64_t mul(uint64_t x, uint64_t y) const { return reduce((u128_local) x * y); }
 uint64_t power(uint64_t x, uint64_t y) {
  uint64_t ret = one;
  for (; y; y >>= 1, x = mul(x, x)) {
   if (y & 1) ret = mul(ret, x);
  }
  return ret;
 }
 uint64_t inv(uint64_t x) { return power(x, mod - 2); }
 uint64_t div_(uint64_t x, uint64_t y) { return mul(x, inv(y)); }
};

inline std::random_device& get_rd() { static std::random_device rd; return rd; }
inline std::mt19937_64& get_rng() { static std::mt19937_64 rng(get_rd()()); return rng; }

inline bool is_prime(uint64_t n) {
 if (n == 1) return false;
 if (n == 2) return true;
 if (n % 2 == 0) return false;
 Montgomery64 mt(n);
 int s = 0;
 uint64_t d = n - 1;
 while (~d & 1) ++s, d >>= 1;
 const int TRY = 5;
 auto& rng = get_rng();
 for (int iter = 0; iter < TRY; ++iter) {
  uint64_t a = rng() % (n - 1) + 1;
  a = mt.power(a, d);
  if (a == mt.one) continue;
  for (int i = 0; i < s && a != mt.sub(0, mt.one); ++i, a = mt.mul(a, a));
  if (a != mt.sub(0, mt.one)) return false;
 }
 return true;
}

inline std::vector<std::pair<uint64_t, int>> factorize_pairs(uint64_t n) {
 if (n == 1) return {};
 std::map<uint64_t, int> mp;
 auto& rng = get_rng();
 auto gen = [&](auto gen, uint64_t n) -> void {
  if (is_prime(n)) {
   ++mp[n];
  } else if (n % 2 == 0) {
   gen(gen, n / 2);
   ++mp[2];
  } else {
   Montgomery64 mt(n);
   while (true) {
    uint64_t a = rng() % n;
    auto f = [&](uint64_t x) -> uint64_t { return mt.add(mt.mul(x, x), a); };
    uint64_t x = rng() % n;
    for (int i = 1;; i <<= 1) {
     uint64_t y = x, q = mt.one;
     for (int j = 0; j < i; ++j) {
      x = f(x);
      q = mt.mul(q, mt.sub(y, x));
     }
     uint64_t g = std::gcd(q, n);
     if (g == 1) continue;
     if (g == n) break;
     gen(gen, g);
     gen(gen, n / g);
     return;
    }
   }
  }
 };
 gen(gen, n);
 return std::vector<std::pair<uint64_t, int>>(mp.begin(), mp.end());
}

inline uint64_t primitive_root(uint64_t p) {
 if (p == 2) return 1;
 auto fact = factorize_pairs(p - 1);
 Montgomery64 mt(p);
 auto& rng = get_rng();
 while (true) {
  uint64_t g = rng() % (p - 1) + 1;
  bool ok = true;
  for (auto [x, y] : fact) {
   if (mt.power(g, (p - 1) / x) == mt.one) { ok = false; break; }
  }
  if (ok) return mt.reduce(g);
 }
}
} // namespace yosupo_336163

struct PrimitiveRoot {
 static vector<u64> run(const vector<u64>& qs) {
  vector<u64> ans;
  ans.reserve(qs.size());
  for (auto p : qs) ans.push_back(yosupo_336163::primitive_root(p));
  return ans;
 }
};
