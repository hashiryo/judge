#pragma once
#include "_common.hpp"
// =============================================================================
// Maspy 流 O(1) mod pow with index calculus precompute
//
// 移植元のアイデア:
//   - LOG table 構築: yosupo-discrete-logarithm-fixed-mod の index_calculus.hpp
//     (purplesyringa 提出 https://judge.yosupo.jp/submission/354636 由来)
//   - 4-stage 構築:
//     Stage 1: 最初の MAGIC4 個 (=25) の小素数の log を Pohlig-Hellman + batch BSGS で
//     Stage 2: その他の小素数を Farey 分解 + 既知 log で再帰計算、smooth な数に伝播
//     Stage 3: spf 線型篩で MAGIC3 までの合成数を埋める
//     Stage 5: MAGIC3+1..MAGIC2 を hos.lyric の漸化式 log(i) = log(-(p%i)) - log(p/i) で
//
//   - 単純 BSGS よりも幅広い p で precompute が高速。
//
// 加えて mod pow 用に:
//   - 順方向 power table T1[i]=g^i, T2[j]=g^{jB} (B = √(p-1)) を追加構築
//   - per-query pow(a, e): log(a) lookup → e * log(a) mod (p-1) → g^x lookup
//
// 制約: p は奇素数, p < 2^32/phi
// =============================================================================
#pragma GCC optimize("O3,unroll-loops")
struct MP {
 using i32 = int32_t;
 static constexpr u32 MAGIC0 = 2000000;
 static constexpr u32 MAGIC1 = 1000000;
 static constexpr u32 MAGIC2 = 1300000;
 static constexpr u32 MAGIC3 = 31624;
 static constexpr u32 MAGIC4 = 25;
 static constexpr u32 MAGIC5 = 100;

 u32 mod;
 u32 mod_minus_1;
 u32 halfp;
 u32 g;
 u32 B;
 std::vector<u32> T1, T2;
 std::vector<u32> dlog_lookup;   // index_calculus と同じ非 centered (size MAGIC2+1 or p)
 std::vector<u32> farey_lookup;

 static u32 inv_mod(u32 a, u32 m) {
  if (a == 0) return 0;
  u32 ori_m = m;
  i32 b = 1, c = 0;
  while (a != 1) {
   u32 d = m / a;
   c -= b * (i32) d;
   m -= a * d;
   std::swap(a, m); std::swap(b, c);
  }
  return b < 0 ? b + ori_m : b;
 }
 static u32 mod_pow(u32 b, i32 e, u32 m) {
  u32 r = 1;
  while (e) { if (e & 1) r = u64(r) * b % m; b = u64(b) * b % m; e >>= 1; }
  return r;
 }
 static u32 crt_combine(u32 a1, u32 m1, u32 a2, u32 m2) {
  return (inv_mod(m1 % m2, m2) * u64(m1) * a2 + inv_mod(m2 % m1, m1) * u64(m2) * a1) % (u64(m1) * m2);
 }

 // Montgomery 風 mulgp_prec (index_calculus.hpp と同じ)
 struct mulgp_prec_t { u32 g_R, neg_inv, p; };
 static mulgp_prec_t mulgp_prec(u32 g, u32 p) {
  u32 ni = p;
  ni *= (2 - ni * p);
  ni *= (2 - ni * p);
  ni *= (2 - ni * p);
  ni *= (2 - ni * p);
  return {u32((u64(g) << 32) % p), u32(-ni), p};
 }
 [[gnu::always_inline]] static inline u32 mulgp(u32 x, const mulgp_prec_t& pr) {
  u64 y = u64(x) * pr.g_R;
  u32 r = u32((y + u64(u32(y) * pr.neg_inv) * pr.p) >> 32);
  return r >= pr.p ? r - pr.p : r;
 }

 struct fast_table {
  u32 cap, mask;
  std::vector<u8> tags;
  std::vector<u32> keys, values;
  fast_table(u32 n) {
   u32 t = (u32) std::ceil(n / 0.8);
   cap = t < 16 ? 16 : 1u << (32 - __builtin_clz(t - 1));
   mask = cap - 1;
   tags.assign(cap, 0); keys.assign(cap, 0); values.assign(cap, u32(-1));
  }
  void insert(u32 key, u32 val) {
   u32 i = key & mask;
   while (tags[i]) i = (i + 1) & mask;
   tags[i] = 1; keys[i] = key; values[i] = val;
  }
  u32 lookup(u32 key) const {
   u32 i = key & mask;
   while (tags[i]) {
    if (keys[i] == key) return values[i];
    i = (i + 1) & mask;
   }
   return u32(-1);
  }
 };

 struct batch_bsgs {
  u32 p, g, order, step_size, invstep;
  mulgp_prec_t invstep_prec;
  fast_table lookup;
  batch_bsgs(u32 p_, u32 g_, u32 order_, u32 ss = 0)
   : p(p_), g(g_), order(order_),
     step_size(ss == 0 ? (u32) std::sqrt(double(order_)) : (ss < 100 ? order_ : ss)),
     lookup(step_size) {
   auto pr = mulgp_prec(g, p);
   for (u32 i = 0, j = 1; i < step_size; ++i, j = mulgp(j, pr)) lookup.insert(j, i);
   invstep = mod_pow(g, order - step_size, p);
   invstep_prec = mulgp_prec(invstep, p);
  }
  u32 solve(u32 h) {
   u32 ret = 0;
   while (lookup.lookup(h) == u32(-1)) {
    ret += step_size;
    h = mulgp(h, invstep_prec);
   }
   return ret + lookup.lookup(h);
  }
 };

 // 原始根探索
 static u32 find_g(u32 p) {
  u32 m = p - 1;
  std::vector<u32> qs;
  u32 t = m;
  for (u32 q = 2; u64(q) * q <= t; ++q) {
   if (t % q == 0) { qs.push_back(q); while (t % q == 0) t /= q; }
  }
  if (t > 1) qs.push_back(t);
  for (u32 g = 2; g < p; ++g) {
   bool ok = true;
   for (u32 q : qs) if (mod_pow(g, m / q, p) == 1) { ok = false; break; }
   if (ok) return g;
  }
  return 0;
 }

 MP() = default;
 MP(u32 m) : mod(m), mod_minus_1(m - 1), halfp(m / 2) {
  g = find_g(m);
  B = u32(std::sqrt(double(mod_minus_1))) + 1;

  // ---- forward power tables T1, T2 ----
  T1.resize(B);
  T1[0] = 1;
  for (u32 i = 1; i < B; ++i) T1[i] = u32(u64(T1[i - 1]) * g % mod);
  u64 gB = u64(T1[B - 1]) * g % mod;
  u32 T2_size = (mod_minus_1 + B - 1) / B + 2;
  T2.resize(T2_size);
  T2[0] = 1;
  for (u32 j = 1; j < T2_size; ++j) T2[j] = u32(u64(T2[j - 1]) * gB % mod);

  // ---- LOG 構築 (index_calculus 流) ----
  if (mod < MAGIC0) {
   dlog_lookup.assign(mod, 0);
   auto pr = mulgp_prec(g, mod);
   u32 i = 1, j = g;
   while (j != 1) { dlog_lookup[j] = i; j = mulgp(j, pr); ++i; }
   return;
  }

  dlog_lookup.assign(MAGIC2 + 1, 0);
  // Farey table
  farey_lookup.assign(MAGIC1, 0);
  {
   auto farey_rec = [&](auto& self, u32 f1, u32 f2, u32 x, u32 y) -> void {
    u32 f3 = f1 + f2;
    u32 l = (((u64) mod * (f3 >> 16) - MAGIC2) * MAGIC1 - 1) / ((u64) mod * (f3 & 0xffff)) + 1;
    u32 r = (((u64) mod * (f3 >> 16) + MAGIC2) * MAGIC1) / ((u64) mod * (f3 & 0xffff));
    l = std::max(l, x); r = std::min(r, y);
    if (x < l) self(self, f1, f3, x, l);
    std::fill(farey_lookup.begin() + l, farey_lookup.begin() + r, f3);
    if (r < y) self(self, f3, f2, r, y);
   };
   const u32 first_x = (u64) MAGIC2 * MAGIC1 / mod;
   const u32 first_y = ((u64) (mod - MAGIC2) * MAGIC1 - 1) / (mod * 2) + 1;
   std::fill(farey_lookup.begin(), farey_lookup.begin() + first_x, 1u);
   farey_rec(farey_rec, 1u, 0x10002u, first_x, first_y);
   std::fill(farey_lookup.begin() + first_y, farey_lookup.begin() + MAGIC1/2, 0x10002u);
   for (u32 i = MAGIC1/2; i < MAGIC1; ++i) {
    farey_lookup[i] = (farey_lookup[MAGIC1 - 1 - i] * 0xffff0001u ^ 0xffff0000u) + 0x10000u;
   }
  }
  // 線型篩 (smallest prime factor)
  std::vector<u32> spf(MAGIC3 + 1, 0), primes;
  primes.reserve(3500);
  for (u32 i = 2; i <= MAGIC3; ++i) {
   if (spf[i] == 0) { spf[i] = i; primes.push_back(i); }
   for (size_t j = 0; primes[j] * i <= MAGIC3; ++j) {
    spf[primes[j] * i] = primes[j];
    if (primes[j] == spf[i]) break;
   }
  }
  // Stage 1: 最初の MAGIC4 個の素数を Pohlig-Hellman で
  std::vector<u32> stage1(primes.begin(), primes.begin() + MAGIC4);
  std::vector<u32> stage1_ans(MAGIC4, 0);
  u32 stage1_mod = 1;
  std::vector<std::pair<u32, u32>> group_size_factorisation;
  {
   u32 gs = mod_minus_1;
   for (u32 pp : primes) {
    if (u64(pp) * pp >= gs) break;
    if (gs % pp == 0) {
     u32 e = 0;
     do { ++e; gs /= pp; } while (gs % pp == 0);
     group_size_factorisation.emplace_back(pp, e);
    }
   }
   if (gs > 1) group_size_factorisation.emplace_back(gs, 1);
  }
  const u32 ginv = inv_mod(g, mod);
  for (auto [pp, e] : group_size_factorisation) {
   u64 nval = 1; for (u32 ii = 0; ii < e; ++ii) nval *= pp;
   u32 ee = u32(mod / nval);
   u32 gg = mod_pow(ginv, ee, mod);
   batch_bsgs ph(mod, mod_pow(g, mod / pp, mod), pp, (u32) std::sqrt(double(pp) * MAGIC4));
   for (u32 i = 0; i < MAGIC4; ++i) {
    u32 h = mod_pow(stage1[i], ee, mod);
    u64 pe = nval; i32 ans = 0; u32 mul_ = 1;
    for (u32 j = 0; j < e; ++j) {
     pe /= pp;
     ans += ph.solve(mod_pow(u64(h) * mod_pow(gg, ans, mod) % mod, u32(pe), mod)) * mul_;
     mul_ *= pp;
    }
    stage1_ans[i] = crt_combine(stage1_ans[i], stage1_mod, ans, u32(nval));
   }
   stage1_mod *= u32(nval);
  }
  for (u32 i = 0; i < MAGIC4; ++i) dlog_lookup[stage1[i]] = stage1_ans[i];

  auto get_frac = [&](u32 h) -> std::pair<i32, u32> {
   u32 bucket = u64(h) * MAGIC1 / mod;
   u32 frac = farey_lookup[bucket];
   u32 den = frac & 0xffff;
   i32 num = (i32) (den * h - (frac >> 16) * mod);
   return {num, den};
  };

  // Stage 2: 残りの小素数 + smooth な数の log
  std::vector<u32> smooth{1};
  smooth.reserve(5000);
  for (u32 pp : primes) {
   if (!dlog_lookup[pp]) {
    u32 h = pp;
    u32 ans = 0;
    while (true) {
     auto [num, den] = get_frac(h);
     u32 anum = (u32) std::abs(num);
     if ((den != 1 && dlog_lookup[den] == 0) || (anum != 1 && dlog_lookup[anum] == 0)) {
      ++ans; h = u64(h) * ginv % mod; continue;
     }
     u32 dlog = (num < 0 ? dlog_lookup[anum] + halfp : dlog_lookup[anum]) + mod_minus_1 - dlog_lookup[den];
     dlog_lookup[pp] = (ans + dlog) % mod_minus_1;
     break;
    }
   }
   if (pp > MAGIC5) continue;
   smooth.erase(std::remove_if(smooth.begin(), smooth.end(),
                               [&](u32 v) { return v > MAGIC2 / pp; }),
                smooth.end());
   u32 nn = (u32) smooth.size();
   for (u32 j = 0; j < nn; ++j) {
    u32 ns = smooth[j] * pp;
    do {
     if (ns <= MAGIC2 / pp) smooth.push_back(ns);
     dlog_lookup[ns] = (dlog_lookup[ns / pp] + dlog_lookup[pp]) % mod_minus_1;
     ns *= pp;
    } while (ns <= MAGIC2);
   }
  }
  // Stage 3: spf で MAGIC3 まで
  for (u32 i = 2; i <= MAGIC3; ++i) {
   if (dlog_lookup[i]) continue;
   dlog_lookup[i] = (dlog_lookup[i / spf[i]] + dlog_lookup[spf[i]]) % mod_minus_1;
  }
  // Stage 5: MAGIC3+1..MAGIC2 を p mod i 漸化式
  for (u32 i = MAGIC3 + 1; i <= MAGIC2; ++i) {
   if (dlog_lookup[i]) continue;
   dlog_lookup[i] = (dlog_lookup[mod % i] + halfp + mod_minus_1 - dlog_lookup[mod / i]) % mod_minus_1;
  }
 }

 // base.cpp 互換 I/F
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
  // log(a) を index_calculus 流 (非 centered) で取得
  u32 logA;
  if (mod < MAGIC0 || a < dlog_lookup.size()) {
   logA = dlog_lookup[a];
  } else {
   u32 bucket = u64(a) * MAGIC1 / mod;
   u32 frac = farey_lookup[bucket];
   u32 den = frac & 0xffff;
   i32 num = (i32) (den * a - (frac >> 16) * mod);
   u32 anum = (u32) std::abs(num);
   u32 r = (num < 0 ? dlog_lookup[anum] + halfp : dlog_lookup[anum]) + mod_minus_1 - dlog_lookup[den];
   if (r >= mod_minus_1) r -= mod_minus_1;
   if (r >= mod_minus_1) r -= mod_minus_1;
   logA = r;
  }
  u64 x = u64(e) * logA % mod_minus_1;
  return u32(u64(T1[x % B]) * T2[x / B] % mod);
 }
};
