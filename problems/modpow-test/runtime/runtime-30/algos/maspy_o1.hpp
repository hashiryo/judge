#pragma once
#include "_common.hpp"
// =============================================================================
// Maspy 流 O(1) mod pow (https://maspypy.com/o1-mod-inv-mod-pow)
// 移植元: https://maspypy.github.io/library/mod/modfast.hpp
//
// 概要:
//   原始根 g を固定し、a^e ≡ g^(e · log_g(a)) (mod p) で計算。
//   - log_g(a) は FRAC table (Farey 分解) + LOG centered table で O(1) lookup
//   - g^x は T1[x mod B] · T2[x / B] (B = √(p-1)) で 2 lookup + 1 mul
//   - per-query: 5 lookup + 2 mul + 2 div ≈ 30-40 cyc
//
// precompute:
//   - 原始根探索: 試し割りで p-1 を素因数分解し、g^((p-1)/q) ≠ 1 の小整数を探す
//   - LOG: 小素数 (i < √p) は BSGS、 大きい i は hos.lyric 漸化式
//          log(i) = log(-(p%i)) - log(p/i)
//   - FRAC: Stern-Brocot 木 BFS (centered.hpp と同じ)
//
// 制約: p < 2^30 (FRAC bucket = a >> 10 を仮定)、p は奇素数
// =============================================================================
struct MP {
 using i32 = int32_t;
 using u16 = unsigned short;
 static constexpr u32 K = 1u << 21;          // LOG table の半幅
 static constexpr u32 FRAC_BUCKETS = 1u << 20;
 static constexpr u32 STERN_LIMIT = 2048;

 u32 mod;
 u32 mod_minus_1;
 u32 g;
 u32 B;                       // ⌈√(p-1)⌉
 vector<u32> T1, T2;          // power tables: g^i, g^(iB)
 vector<u32> LOG;             // [-K, K] dlog, centered
 vector<pair<u16, u16>> FRAC; // bucket → (num, den)

 static u32 pow_mod(u64 a, u64 e, u32 p) {
  u64 r = 1;
  while (e) { if (e & 1) r = r * a % p; a = a * a % p; e >>= 1; }
  return u32(r);
 }

 static u32 find_g(u32 p) {
  u32 m = p - 1;
  vector<u32> qs;
  u32 t = m;
  for (u32 q = 2; u64(q) * q <= t; ++q) {
   if (t % q == 0) { qs.push_back(q); while (t % q == 0) t /= q; }
  }
  if (t > 1) qs.push_back(t);
  for (u32 g = 2; g < p; ++g) {
   bool ok = true;
   for (u32 q : qs) if (pow_mod(g, m / q, p) == 1) { ok = false; break; }
   if (ok) return g;
  }
  return 0;
 }

 MP() = default;
 MP(u32 m) : mod(m), mod_minus_1(m - 1) {
  g = find_g(m);
  B = u32(std::sqrt(double(mod_minus_1))) + 1;

  // ---- T1, T2 構築 ----
  T1.resize(B);
  T1[0] = 1;
  for (u32 i = 1; i < B; ++i) T1[i] = u32(u64(T1[i - 1]) * g % mod);
  u64 gB = u64(T1[B - 1]) * g % mod;
  u32 T2_size = (mod_minus_1 + B - 1) / B + 2;
  T2.resize(T2_size);
  T2[0] = 1;
  for (u32 j = 1; j < T2_size; ++j) T2[j] = u32(u64(T2[j - 1]) * gB % mod);

  // ---- BSGS dlog 用 sorted baby table ----
  vector<pair<u32, u32>> baby(B);
  for (u32 i = 0; i < B; ++i) baby[i] = {T1[i], i};
  std::sort(baby.begin(), baby.end());
  u32 inv_gB = pow_mod(gB, mod - 2, mod);

  auto bsgs = [&](u64 x) -> u32 {
   u64 cur = x;
   u32 max_a = u32((u64(mod_minus_1) + B - 1) / B) + 1;
   for (u32 a = 0; a <= max_a; ++a) {
    auto it = std::lower_bound(baby.begin(), baby.end(), std::make_pair(u32(cur), u32(0)));
    if (it != baby.end() && it->first == cur) {
     u64 r = u64(a) * B + it->second;
     return u32(r % mod_minus_1);
    }
    cur = cur * inv_gB % mod;
   }
   return u32(-1);
  };

  // ---- LOG table 構築 ----
  LOG.assign(2 * K + 1, 0);
  LOG[K + 1] = 0;
  // log(-1) = (p-1)/2
  // (これは LOG[K - 1] = (p-1)/2、つまり後段の負側 fill で確定)

  u32 sqrt_p = u32(std::sqrt(double(mod))) + 1;

  // 小さい i ∈ [2, min(K, sqrt_p)) は素数のみ BSGS、合成数は線型篩で
  vector<u32> spf(std::min(K, sqrt_p) + 1, 0);
  vector<u32> primes;
  u32 sp_lim = std::min(K, sqrt_p);
  for (u32 i = 2; i <= sp_lim; ++i) {
   if (spf[i] == 0) { spf[i] = i; primes.push_back(i); }
   for (u32 q : primes) {
    if (q > spf[i] || u64(q) * i > sp_lim) break;
    spf[q * i] = q;
   }
  }
  for (u32 i = 2; i <= sp_lim; ++i) {
   if (spf[i] == i) {
    LOG[K + i] = bsgs(i);
   } else {
    u32 a = LOG[K + spf[i]];
    u32 b = LOG[K + i / spf[i]];
    LOG[K + i] = u32((u64(a) + b) % mod_minus_1);
   }
  }

  // 大きい i ∈ [sqrt_p, K]: hos.lyric 漸化式
  // log(i) = log(-(p%i)) - log(p/i) = (p-1)/2 + log(p%i) - log(p/i)
  for (u32 i = sqrt_p; i <= K; ++i) {
   u32 r = mod % i;
   u32 q = mod / i;
   u32 logR = LOG[K + r];          // log(p%i)
   u32 logQ = LOG[K + q];          // log(p/i)
   // (p-1)/2 + logR - logQ mod (p-1)
   u64 v = u64(mod_minus_1 / 2) + logR;
   if (v >= mod_minus_1) v -= mod_minus_1;
   u32 vv = u32(v) >= logQ ? u32(v) - logQ : u32(v) + mod_minus_1 - logQ;
   LOG[K + i] = vv;
  }

  // 負側: LOG[K - i] = (p-1)/2 + LOG[K + i] mod (p-1)
  for (u32 i = 1; i <= K; ++i) {
   u64 v = u64(LOG[K + i]) + mod_minus_1 / 2;
   if (v >= mod_minus_1) v -= mod_minus_1;
   LOG[K - i] = u32(v);
  }

  // ---- FRAC table (centered.hpp と同じ Stern-Brocot BFS) ----
  FRAC.assign(FRAC_BUCKETS + 1, {0, 0});
  std::vector<std::tuple<u32, u32, u32, u32>> stk;
  stk.emplace_back(0, 1, 1, 1);
  while (!stk.empty()) {
   auto [aa, bb, cc, dd] = stk.back(); stk.pop_back();
   if (bb + dd < STERN_LIMIT) {
    stk.emplace_back(aa + cc, bb + dd, cc, dd);
    stk.emplace_back(aa, bb, aa + cc, bb + dd);
    continue;
   }
   u32 s = u32(u64(aa) * mod / (1024u * bb));
   u32 t = u32(u64(cc) * mod / (1024u * dd));
   if (s <= FRAC_BUCKETS) FRAC[s] = {u16(aa), u16(bb)};
   if (t <= FRAC_BUCKETS) FRAC[t] = {u16(cc), u16(dd)};
   const u32 a_min = std::min(aa, cc), b_min = std::min(bb, dd);
   for (u32 i = s + 1; i < t && i <= FRAC_BUCKETS; ++i) {
    FRAC[i] = {u16(a_min), u16(b_min)};
   }
  }
 }

 // base.cpp 互換の I/F
 inline u32 set(u32 n) const { return n; }
 inline u32 get(u32 n) const { return n; }
 inline u32 norm(u32 n) const { return n; }
 inline u32 mul(u32 l, u32 r) const { return u32(u64(l) * r % mod); }
 inline u32 plus(u32 l, u32 r) const { return l += r, l < mod ? l : l - mod; }
 inline u32 diff(u32 l, u32 r) const { return l -= r, l >> 31 ? l + mod : l; }

 inline u32 pow(u32 a, u32 e) const {
  if (a == 0) return e == 0 ? 1 : 0;
  if (e == 0) return 1;
  if (a == 1) return 1;
  // log(a) を FRAC + LOG で取得
  const auto [num, den] = FRAC[a >> 10];
  const u32 t = a * den - u32(num) * mod;  // u32 wrap で signed |t| ≤ K
  const u32 idx_num = K + t;               // wrap で正しい centered index
  const u32 idx_den = K + den;
  const u32 La = LOG[idx_num];
  const u32 Lb = LOG[idx_den];
  const u32 logA = La >= Lb ? La - Lb : La + mod_minus_1 - Lb;
  // x = e · logA mod (p-1)
  const u64 x = u64(e) * logA % mod_minus_1;
  // g^x = T1[x mod B] · T2[x / B]
  return u32(u64(T1[x % B]) * T2[x / B] % mod);
 }
};
