#pragma once
#include "_common.hpp"
// =============================================================================
// Maspy 流 index calculus: 「rational approximation + factor base」で per-query O(1)。
//
// 移植元: yosupo discrete_logarithm_fixed_mod の最速提出 (purplesyringa, 2026-04)
//   https://judge.yosupo.jp/submission/354636
// アルゴリズム解説:
//   https://maspypy.com/o1-mod-inv-mod-pow
//
// アイデア:
//   x ∈ Z/p^* について、Stern-Brocot/Farey で x ≡ a/b mod p, |a|, b ≤ ~√p に表現。
//   小素数 (factor base) の log を事前計算 → a, b を small primes に factor → log 加算。
//
// per query: ~O(log p) (Farey 1 lookup + factor 加算)
// precompute: 100-300 ms (factor base log 計算 + Farey 表構築)
//
// 実装定数 (purplesyringa 提出と同じ):
//   MAGIC0 = 2000000   // p < これなら全 log を直接 precompute
//   MAGIC1 = 1000000   // Farey buckets
//   MAGIC2 = 1300000   // Farey で得られる |a| / b の上限 ≈
//   MAGIC3 = 31624     // factorisation precompute 上限
//   MAGIC4 = 25        // batch dlog する小素数の個数
//   MAGIC5 = 100       // smooth net 拡張閾値
//
// 注: 元提出は blazingio という超高速 I/O ライブラリも使用していたが、ここでは
//     アルゴリズム部分のみを移植 (I/O は base.cpp の標準 cin/cout)。
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
struct DLog {
 using i32 = int32_t;
 using fast_u32 = u32;

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

 // Montgomery 風 mulgp_prec
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

 // 開放アドレス hash table (BSGS 用)
 struct fast_table {
  u32 cap, mask;
  std::vector<u8> tags;
  std::vector<u32> keys, values;
  fast_table(u32 n) {
   u32 t = (u32) std::ceil(n / 0.8);
   cap = t < 16 ? 16 : 1u << (32 - __builtin_clz(t - 1));
   mask = cap - 1;
   tags.assign(cap, 0);
   keys.assign(cap, 0);
   values.assign(cap, u32(-1));
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

 // BSGS subset (Pohlig-Hellman 内で使う)
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

 // 主クラス: fixed-mod 版 (precompute 1 回 + per-query O(1))
 // 実装は purplesyringa 提出 (https://judge.yosupo.jp/submission/354636) を移植。
 static constexpr u32 MAGIC0 = 2000000;
 static constexpr u32 MAGIC1 = 1000000;
 static constexpr u32 MAGIC2 = 1300000;
 static constexpr u32 MAGIC3 = 31624;
 static constexpr u32 MAGIC4 = 25;
 static constexpr u32 MAGIC5 = 100;

 static std::vector<u32> run(u32 p, u32 g, const std::vector<u32>& qs) {
  const u32 pm1 = p - 1;
  const u32 halfp = p / 2;
  std::vector<u32> dlog_lookup;
  std::vector<u32> farey_lookup;

  // 全部 precompute (p が小さい場合)
  if (p < MAGIC0) {
   dlog_lookup.assign(p, 0);
   auto pr = mulgp_prec(g, p);
   u32 i = 1, j = g;
   while (j != 1) { dlog_lookup[j] = i; j = mulgp(j, pr); ++i; }
  } else {
   dlog_lookup.assign(MAGIC2 + 1, 0);
   // Farey fractions
   farey_lookup.assign(MAGIC1, 0);
   {
    auto farey_rec = [&](auto& self, u32 f1, u32 f2, u32 x, u32 y) -> void {
     u32 f3 = f1 + f2;
     u32 l = (((u64) p * (f3 >> 16) - MAGIC2) * MAGIC1 - 1) / ((u64) p * (f3 & 0xffff)) + 1;
     u32 r = (((u64) p * (f3 >> 16) + MAGIC2) * MAGIC1) / ((u64) p * (f3 & 0xffff));
     l = std::max(l, x); r = std::min(r, y);
     if (x < l) self(self, f1, f3, x, l);
     std::fill(farey_lookup.begin() + l, farey_lookup.begin() + r, f3);
     if (r < y) self(self, f3, f2, r, y);
    };
    const u32 first_x = (u64) MAGIC2 * MAGIC1 / p;
    const u32 first_y = ((u64) (p - MAGIC2) * MAGIC1 - 1) / (p * 2) + 1;
    std::fill(farey_lookup.begin(), farey_lookup.begin() + first_x, 1u);
    farey_rec(farey_rec, 1u, 0x10002u, first_x, first_y);
    std::fill(farey_lookup.begin() + first_y, farey_lookup.begin() + MAGIC1/2, 0x10002u);
    for (u32 i = MAGIC1/2; i < MAGIC1; ++i) {
     farey_lookup[i] = (farey_lookup[MAGIC1 - 1 - i] * 0xffff0001u ^ 0xffff0000u) + 0x10000u;
    }
   }
   // factorisation precompute up to MAGIC3
   std::vector<u32> spf(MAGIC3 + 1, 0), primes;
   primes.reserve(3500);
   for (u32 i = 2; i <= MAGIC3; ++i) {
    if (spf[i] == 0) { spf[i] = i; primes.push_back(i); }
    for (size_t j = 0; primes[j] * i <= MAGIC3; ++j) {
     spf[primes[j] * i] = primes[j];
     if (primes[j] == spf[i]) break;
    }
   }
   // Stage 1: 最初の MAGIC4 個の素数の log を Pohlig-Hellman で求める
   std::vector<u32> stage1(primes.begin(), primes.begin() + MAGIC4);
   std::vector<u32> stage1_ans(MAGIC4, 0);
   u32 stage1_mod = 1;
   std::vector<std::pair<u32, u32>> group_size_factorisation;
   {
    u32 gs = pm1;
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
   const u32 ginv = inv_mod(g, p);
   for (auto [pp, e] : group_size_factorisation) {
    u32 n = mod_pow(pp, e, p); (void) n;  // not used directly
    u64 nval = 1; for (u32 ii = 0; ii < e; ++ii) nval *= pp;
    u32 ee = u32(p / nval);
    u32 gg = mod_pow(ginv, ee, p);
    batch_bsgs ph(p, mod_pow(g, p / pp, p), pp, (u32) std::sqrt(double(pp) * MAGIC4));
    for (u32 i = 0; i < MAGIC4; ++i) {
     u32 h = mod_pow(stage1[i], ee, p);
     u64 pe = nval; i32 ans = 0; u32 mul_ = 1;
     for (u32 j = 0; j < e; ++j) {
      pe /= pp;
      ans += ph.solve(mod_pow(u64(h) * mod_pow(gg, ans, p) % p, u32(pe), p)) * mul_;
      mul_ *= pp;
     }
     stage1_ans[i] = crt_combine(stage1_ans[i], stage1_mod, ans, u32(nval));
    }
    stage1_mod *= u32(nval);
   }
   for (u32 i = 0; i < MAGIC4; ++i) dlog_lookup[stage1[i]] = stage1_ans[i];

   // get_frac
   auto get_frac = [&](u32 h) -> std::pair<i32, u32> {
    u32 bucket = u64(h) * MAGIC1 / p;
    u32 frac = farey_lookup[bucket];
    u32 den = frac & 0xffff;
    i32 num = (i32) (den * h - (frac >> 16) * p);
    return {num, den};
   };

   // Stage 2: その他の小素数 + smooth number の log
   std::vector<u32> smooth{1};
   smooth.reserve(5000);
   for (u32 pp : primes) {
    if (!dlog_lookup[pp]) {
     // 2a: pp の log を rational 表現で計算
     u32 h = pp;
     u32 ans = 0;
     while (true) {
      auto [num, den] = get_frac(h);
      u32 anum = (u32) std::abs(num);
      if ((den != 1 && dlog_lookup[den] == 0) || (anum != 1 && dlog_lookup[anum] == 0)) {
       ++ans; h = u64(h) * ginv % p; continue;
      }
      u32 dlog = (num < 0 ? dlog_lookup[anum] + halfp : dlog_lookup[anum]) + pm1 - dlog_lookup[den];
      dlog_lookup[pp] = (ans + dlog) % pm1;
      break;
     }
    }
    if (pp > MAGIC5) continue;
    // 2b: smooth な数の log を波及計算
    smooth.erase(std::remove_if(smooth.begin(), smooth.end(),
                                [&](u32 v) { return v > MAGIC2 / pp; }),
                 smooth.end());
    u32 nn = (u32) smooth.size();
    for (u32 j = 0; j < nn; ++j) {
     u32 ns = smooth[j] * pp;
     do {
      if (ns <= MAGIC2 / pp) smooth.push_back(ns);
      dlog_lookup[ns] = (dlog_lookup[ns / pp] + dlog_lookup[pp]) % pm1;
      ns *= pp;
     } while (ns <= MAGIC2);
    }
   }
   // Stage 3: spf を使って MAGIC3 までの全数の log
   for (u32 i = 2; i <= MAGIC3; ++i) {
    if (dlog_lookup[i]) continue;
    dlog_lookup[i] = (dlog_lookup[i / spf[i]] + dlog_lookup[spf[i]]) % pm1;
   }
   // Stage 5: MAGIC3+1 から MAGIC2 までを p mod i トリックで
   for (u32 i = MAGIC3 + 1; i <= MAGIC2; ++i) {
    if (dlog_lookup[i]) continue;
    dlog_lookup[i] = (dlog_lookup[p % i] + halfp + pm1 - dlog_lookup[p / i]) % pm1;
   }
  }

  // Per-query lookup
  std::vector<u32> ans;
  ans.reserve(qs.size());
  if (p < MAGIC0) {
   for (u32 a : qs) {
    if (a == 0) ans.push_back(u32(-1));
    else if (a == 1) ans.push_back(0);
    else ans.push_back(dlog_lookup[a]);
   }
  } else {
   for (u32 a : qs) {
    if (a == 0) { ans.push_back(u32(-1)); continue; }
    if (a < dlog_lookup.size()) { ans.push_back(dlog_lookup[a]); continue; }
    u32 bucket = u64(a) * MAGIC1 / p;
    u32 frac = farey_lookup[bucket];
    u32 den = frac & 0xffff;
    i32 num = (i32) (den * a - (frac >> 16) * p);
    u32 anum = (u32) std::abs(num);
    u32 r = (num < 0 ? dlog_lookup[anum] + halfp : dlog_lookup[anum]) + pm1 - dlog_lookup[den];
    if (r >= pm1) r -= pm1;
    if (r >= pm1) r -= pm1;
    ans.push_back(r);
   }
  }
  return ans;
 }
};
