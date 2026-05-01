#pragma once
// =============================================================================
// 試作: yosupo_fastest.hpp (Pollard-Rho with Brent + GCD batching + 2-track)
//   をベースに、以下を差し替えた hybrid:
//     - 素数判定: 7-base MR → FJ64_262K (yosupo-primality-test の yosupo_fastest)
//                hash で base 1 個引いて 2-base MR (table 512KB)
//     - GCD: std::gcd (Euclid) → mylib の binary_gcd
//   どちらも Pollard-Rho の内ループには直接効かないが、batch 終端で呼ぶ
//   gcd と Is_Prime の定数を削る試み。
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include <concepts>
#include "mylib/number_theory/binary_gcd.hpp"

// FJ64_262K の prime_test を流用 (primality-test 側の yosupo_fastest.hpp と同一)。
// 名前空間を分けて持ち込む。
namespace fj64_for_factorize {
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;
using u128_local = __uint128_t;

namespace bit_local {
inline int low(uint64_t x) { return x == 0 ? -1 : __builtin_ctzll(x); }
}

// bases table 262144 (u16). primality-test/algos/yosupo_fastest.hpp と同一データ。
// ファイル肥大化を避けるため、共有テーブルは prime test のヘッダから include する。
} // namespace fj64_for_factorize

// FJ64_262K table + prime_test を取り込む (primality-test の fastest を再利用)。
// algos/ 内の相対 include はテストハーネスでも -I が通るので問題なし。
#include "../../yosupo-primality-test/algos/yosupo_fastest.hpp"
// 上で algos/_common.hpp が再 include されるが #pragma once で防がれる。
// `Primality` struct も定義されてしまうのでこの hpp 単体ではコンパイル不能。
// → algos/yosupo_fastest.hpp 内で Primality を匿名化したいが侵襲的なので、
//   ここでは prime_test を別 namespace から呼ぶだけにする。

namespace yosupo_factorize_mixed {
using std::uint32_t;
using std::uint64_t;

template<class U0, class U1>
struct Montgomery {
 constexpr static unsigned B0 = sizeof(U0) * 8U;
 U0 n, nr, rs, np;
 constexpr Montgomery(const U0& Mod) { SetMod(Mod); }
 constexpr void SetMod(const U0& Mod) {
  n = nr = Mod, rs = -static_cast<U1>(n) % n;
  for (uint32_t i = 0; i < std::bit_width<unsigned>(B0) - 1; ++i) nr *= 2 - n * nr;
  np = Reduce(static_cast<U0>(1), rs);
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
 constexpr U0 mul(const U0& x, const U0& y) const noexcept { return Reduce(x, y); }
 constexpr U0 mul_add(const U0& x, const U0& y, const U0& z) const noexcept { return Reduce(x, y, z); }
};

template<bool sorted>
inline std::vector<std::pair<uint64_t, uint32_t>> Factorize(uint64_t n) {
 std::vector<std::pair<uint64_t, uint32_t>> ans;
 if (n % 2 == 0) {
  uint32_t z = __builtin_ctzll(n);
  ans.push_back({2ULL, z}), n >>= z;
 }
 auto upd = [&](const uint64_t& x) {
  for (auto& [p, c] : ans) if (x == p) { ++c; return; }
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
     // gcd を std::gcd → binary_gcd に差し替え
     const uint64_t q3 = Mod.mul(q1, q2), g3 = binary_gcd<uint64_t>(nn, q3);
     if (g3 == 1) continue;
     if (g3 != nn) { ans2 = g3; return; }
     const uint64_t g1 = binary_gcd<uint64_t>(nn, q1);
     const uint64_t g2 = binary_gcd<uint64_t>(nn, q2);
     const uint64_t C = g1 != 1 ? C1 : C2;
     const uint64_t x = g1 != 1 ? x1 : x2;
     uint64_t z = g1 != 1 ? y1 : y2;
     uint64_t g = g1 != 1 ? g1 : g2;
     if (g == nn) {
      do {
       z = Mod.mul_add(z, z, C);
       g = binary_gcd<uint64_t>(nn, x - z);
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
  // 素数判定を FJ64_262K に差し替え
  if (yosupo_fj64::nt::prime_test(nn)) return upd(nn);
  uint64_t d = Pollard_Rho(nn);
  self(self, d), self(self, nn / d);
 };
 if (n > 1) DFS(DFS, n);
 if constexpr (sorted) std::sort(ans.begin(), ans.end());
 return ans;
}
} // namespace yosupo_factorize_mixed

// primality-test/algos/yosupo_fastest.hpp が定義した Primality を捨てて、
// この hpp 用の Factorize を定義する。
#undef ALGO_HPP
struct Factorize {
 static vector<vector<u64>> run(const vector<u64>& qs) {
  vector<vector<u64>> ans;
  ans.reserve(qs.size());
  for (auto x : qs) {
   vector<u64> fs;
   if (x > 1) {
    auto f = yosupo_factorize_mixed::Factorize<true>(x);
    for (auto [p, c] : f) {
     for (uint32_t k = 0; k < c; ++k) fs.push_back(p);
    }
   }
   ans.push_back(std::move(fs));
  }
  return ans;
 }
};
