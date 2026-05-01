#pragma once
// =============================================================================
// Source: yosupo "Factorize" 提出 317402.
//   https://judge.yosupo.jp/submission/317402
//   方式: Brent's Pollard-Rho + Montgomery64 で 64-bit 整数を素因数分解。
//     - 7-base Miller-Rabin (Jaeschke set) で素数判定
//     - 合成数なら Pollard-Rho で 1 つの非自明因子を抽出 → 再帰
//     - Brent ロー検出 + GCD バッチング (M=600 で gcd 呼び出しを 1/M に削減)
//     - 2 トラック並列 (C1=1, C2=2 の独立 Rho を同時に走らせる) で運用
//   Montgomery struct は templated (U0=u64, U1=u128) で `mul_add(x,y,z)` =
//   `Reduce(x*y) + z` を 1 命令でできる API になっている (gcd 後の累積に便利)。
// 抽出方針:
//   - 元提出の I/O (rd/wr) は base.cpp に移譲
//   - Montgomery + Is_Prime + Factorize<sorted=true> はそのままコピー
//   - Factorize::run wrapper で各クエリを vector<u64> に展開
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include <concepts>

namespace yosupo_factorize_317402 {
using std::uint32_t;
using std::uint64_t;

template<class U0, class U1>
struct Montgomery {
 constexpr static unsigned B0 = sizeof(U0) * 8U;
 U0 n, nr, rs, np;

 constexpr Montgomery(const U0& Mod) { SetMod(Mod); }

 constexpr U0 GetMod() const noexcept { return n; }
 constexpr void SetMod(const U0& Mod) {
  assert(Mod >= 2), assert(Mod % 2 == 1);
  assert((Mod >> (B0 - 2)) == 0);
  n = nr = Mod, rs = -static_cast<U1>(n) % n;
  // libc++ には std::__lg が無いので bit_width 経由 (= floor log2)
  for (uint32_t i = 0; i < std::bit_width<unsigned>(B0) - 1; ++i) nr *= 2 - n * nr;
  np = Reduce(static_cast<U0>(1), rs);
 }
 constexpr U0 Reduce(const U0& x) const noexcept {
  const U0 q = x * nr;
  const U0 m = (static_cast<U1>(q) * n) >> B0;
  return n - m;
 }
 constexpr U0 Reduce(const U0& x, const U0& y) const noexcept {
  const U1 t = static_cast<U1>(x) * y;
  const U0 c = (U0) t, d = (U0) (t >> B0);
  const U0 q = c * nr;
  const U0 m = (static_cast<U1>(q) * n) >> B0;
  return d + n - m;
 }
 constexpr U0 Reduce(const U0& x, const U0& y, const U0& z) const noexcept {
  const U1 t = static_cast<U1>(x) * y;
  const U0 c = (U0) t, d = (U0) (t >> B0);
  const U0 q = c * nr;
  const U0 m = (static_cast<U1>(q) * n) >> B0;
  return z + d + n - m;
 }
 constexpr U0 val(const U0& x) const noexcept {
  const uint64_t t = Reduce(x);
  return (t == n) ? static_cast<U0>(0) : t;
 }
 constexpr U0 zero() const noexcept { return static_cast<U0>(0); }
 constexpr U0 one() const noexcept { return np; }
 constexpr U0 raw(const U0& x) const noexcept { return Reduce(x, rs); }
 template<class U> requires std::unsigned_integral<U>
 constexpr U0 trans(const U& x) const noexcept {
  if (__builtin_expect(x < n, 1)) return raw(x);
  return Reduce(x % n, rs);
 }
 constexpr U0 add(const U0& x, const U0& y) const noexcept {
  return (x + y >= 2 * n) ? (x + y - 2 * n) : (x + y);
 }
 constexpr U0 sub(const U0& x, const U0& y) const noexcept {
  return (x < y) ? (x - y + 2 * n) : (x - y);
 }
 constexpr U0 mul(const U0& x, const U0& y) const noexcept { return Reduce(x, y); }
 constexpr U0 mul_add(const U0& x, const U0& y, const U0& z) const noexcept { return Reduce(x, y, z); }
 constexpr bool same(const U0& x, const U0& y) const noexcept {
  const U0 dif = x - y;
  return (dif == 0) || (dif == n) || (dif == -n);
 }
};

constexpr bool Is_Prime(uint64_t x) noexcept {
 if (x <= 1) return false;
 if (x % 2 == 0) return x == 2;
 constexpr std::array<uint64_t, 11> Base{2, 3, 5, 7, 2, 325, 9375, 28178, 450775, 9780504, 1795265022};
 const uint32_t s = __builtin_ctzll(x - 1);
 const uint64_t d = (x - 1) >> s;
 const int q = 63 ^ __builtin_clzll(d);
 const Montgomery<uint64_t, __uint128_t> Mod(x);
 const int l = (x >> 32) ? 4 : 0;
 const int r = (x >> 32) ? 11 : 4;
 for (int _ = l; _ < r; ++_) {
  uint64_t base = Base[_];
  if (base % x == 0) continue;
  base = Mod.trans(base);
  uint64_t a = base;
  for (int i = q - 1; ~i; --i) {
   a = Mod.mul(a, a);
   if ((d >> i) & 1) a = Mod.mul(a, base);
  }
  if (Mod.same(a, Mod.one())) continue;
  for (uint32_t t = 1; t < s && !Mod.same(a, x - Mod.one()); ++t) a = Mod.mul(a, a);
  if (!Mod.same(a, x - Mod.one())) return false;
 }
 return true;
}

template<bool sorted>
inline std::vector<std::pair<uint64_t, uint32_t>> Factorize(uint64_t n) {
 std::vector<std::pair<uint64_t, uint32_t>> ans;
 if (n % 2 == 0) {
  uint32_t z = __builtin_ctzll(n);
  ans.push_back({2ULL, z}), n >>= z;
 }
 auto upd = [&](const uint64_t& x) {
  for (auto& [p, c] : ans) {
   if (x == p) { ++c; return; }
  }
  ans.push_back({x, 1});
 };
 auto Pollard_Rho = [&](const uint64_t& nn) -> uint64_t {
  if (nn % 2 == 0) return 2ULL;
  const Montgomery<uint64_t, __uint128_t> Mod(nn);
  const uint64_t C1 = 1, C2 = 2, M = 600;
  uint64_t Z1 = 1, Z2 = 2, ans2 = 0;
  auto find = [&]() {
   uint64_t z1 = Z1, z2 = Z2;
   for (uint64_t k = M;; k *= 2) {
    const uint64_t x1 = z1 + nn, x2 = z2 + nn;
    for (uint64_t j = 0; j < k; j += M) {
     const uint64_t y1 = z1, y2 = z2;
     uint64_t q1 = 1, q2 = 2;
     z1 = Mod.mul_add(z1, z1, C1), z2 = Mod.mul_add(z2, z2, C2);
     for (uint64_t i = 0; i < M; ++i) {
      uint64_t t1 = x1 - z1, t2 = x2 - z2;
      z1 = Mod.mul_add(z1, z1, C1), z2 = Mod.mul_add(z2, z2, C2);
      q1 = Mod.mul(q1, t1), q2 = Mod.mul(q2, t2);
     }
     q1 = Mod.mul(q1, x1 - z1), q2 = Mod.mul(q2, x2 - z2);
     const uint64_t q3 = Mod.mul(q1, q2), g3 = std::gcd(nn, q3);
     if (g3 == 1) continue;
     if (g3 != nn) { ans2 = g3; return; }
     const uint64_t g1 = std::gcd(nn, q1);
     const uint64_t g2 = std::gcd(nn, q2);
     const uint64_t C = g1 != 1 ? C1 : C2;
     const uint64_t x = g1 != 1 ? x1 : x2;
     uint64_t z = g1 != 1 ? y1 : y2;
     uint64_t g = g1 != 1 ? g1 : g2;
     if (g == nn) {
      do {
       z = Mod.mul_add(z, z, C);
       g = std::gcd(nn, x - z);
      } while (g == 1);
     }
     if (g != nn) { ans2 = g; return; }
     Z1 += 2, Z2 += 2;
     return;
    }
   }
  };
  do { find(); } while (!ans2);
  return ans2;
 };
 auto DFS = [&](auto&& self, const uint64_t& nn) -> void {
  if (Is_Prime(nn)) return upd(nn);
  uint64_t d = Pollard_Rho(nn);
  self(self, d), self(self, nn / d);
 };
 if (n > 1) DFS(DFS, n);
 if constexpr (sorted) std::sort(ans.begin(), ans.end());
 return ans;
}
} // namespace yosupo_factorize_317402

struct Factorize {
 static vector<vector<u64>> run(const vector<u64>& qs) {
  vector<vector<u64>> ans;
  ans.reserve(qs.size());
  for (auto x : qs) {
   vector<u64> fs;
   if (x > 1) {
    auto f = yosupo_factorize_317402::Factorize<true>(x);
    for (auto [p, c] : f) {
     for (uint32_t k = 0; k < c; ++k) fs.push_back(p);
    }
   }
   ans.push_back(std::move(fs));
  }
  return ans;
 }
};
