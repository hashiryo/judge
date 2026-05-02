#pragma once
// 素朴 reference: BSGS (baby-step giant-step)。位数 2^64 - 1 で動く。
//   - baby:  m = ceil(sqrt(N))。N = 2^64 - 1 ≈ 1.8e19、m ≈ 4.3e9 → メモリ不足
// なので**素朴 BSGS は無理**。Pohlig-Hellman で部分群分解しないと正解できない。
//
// このファイルは「正しいけど遅い」reference として、Pohlig-Hellman で
// 各素因数部分群の BSGS を実行し、CRT で合成する版。
//
// 2^64 - 1 = 3 × 5 × 17 × 257 × 641 × 65537 × 6700417
// 最大素因数 6700417 → 部分群 BSGS で sqrt(6.7M) ≈ 2588 baby step。これなら可。
#include "_common.hpp"

namespace gf2_64_ref {
inline std::pair<u64, u64> clmul_loop(u64 a, u64 b) {
 u64 lo = 0, hi = 0;
 for (int i = 0; i < 64; ++i) {
  if ((b >> i) & 1) { lo ^= a << i; if (i) hi ^= a >> (64 - i); }
 }
 return {lo, hi};
}
inline u64 reduce_naive(u64 lo, u64 hi) {
 for (int i = 63; i >= 0; --i) {
  if ((hi >> i) & 1) {
   hi ^= u64(1) << i;
   lo ^= IRRED_LOW << i;
   if (i > 0) hi ^= IRRED_LOW >> (64 - i);
  }
 }
 return lo;
}
inline u64 mul(u64 a, u64 b) { auto [lo, hi] = clmul_loop(a, b); return reduce_naive(lo, hi); }
inline u64 pow(u64 a, u64 e) {
 u64 r = 1;
 while (e) { if (e & 1) r = mul(r, a); a = mul(a, a); e >>= 1; }
 return r;
}

// 2^64 - 1 の素因数
constexpr std::array<u64, 7> ORDER_PRIMES = {3, 5, 17, 257, 641, 65537, 6700417};

// 部分群 (位数 q) 内での BSGS。base, target は既にその部分群に正規化されている前提。
// log_base(target) を [0, q) で返す。
inline u64 bsgs_subgroup(u64 base, u64 target, u64 q) {
 // m = ceil(sqrt(q))
 u64 m = 1;
 while (m * m < q) ++m;
 std::unordered_map<u64, u64> table;
 table.reserve(m * 2);
 u64 cur = 1;
 for (u64 j = 0; j < m; ++j) {
  // 既に存在する j があるなら更新しない (最小の j を使う)
  table.try_emplace(cur, j);
  cur = mul(cur, base);
 }
 // giant step: base^(-m) = base^(q - m) (既に q が位数なので)
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
 return ~u64(0); // 見つからない (起こらないはず)
}

inline u64 log_g(u64 x) {
 // Pohlig-Hellman: 各素因数 q_i ごとに log mod q_i を BSGS で求めて CRT。
 // (2^64 - 1 は square-free なので素因数 1 個ずつで OK)
 u64 N = GROUP_ORDER;
 u64 g = LOG_GENERATOR;
 u64 result = 0; // 最終的に mod N で合成
 u64 mod = 1;    // CRT の現在の moduli 積
 for (u64 q : ORDER_PRIMES) {
  u64 g_sub = pow(g, N / q);   // 部分群 (位数 q) の生成元
  u64 x_sub = pow(x, N / q);   // x をその部分群に射影
  u64 r = bsgs_subgroup(g_sub, x_sub, q); // log mod q
  // CRT: result mod mod を、新しい合成: result' mod (mod * q) で更新
  // result + mod * t ≡ r (mod q) を解く。t = (r - result) * mod^(-1) mod q
  // ここでは mod は q たちの積、各 q は互いに素なので逆元が一意。
  u64 inv_mod_q = 1;
  // mod^(-1) mod q を Fermat (q は素数) で
  u64 mm = mod % q, e = q - 2;
  while (e) { if (e & 1) inv_mod_q = (__uint128_t) inv_mod_q * mm % q; mm = (__uint128_t) mm * mm % q; e >>= 1; }
  u64 diff = (r >= result % q) ? (r - result % q) : (q - (result % q - r));
  u64 t = (__uint128_t) diff * inv_mod_q % q;
  result += mod * t;
  mod *= q;
 }
 return result;
}
} // namespace gf2_64_ref

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& xs) {
  vector<u64> ans(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) ans[i] = gf2_64_ref::log_g(xs[i]);
  return ans;
 }
};
