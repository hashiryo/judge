#pragma once
#include "_common.hpp"
// =============================================================================
// maspy_o1_plantard.hpp の inv_lookup 構築を Plantard 乗算ベースの offline batch
// で置換。
//
// 鍵: Montgomery 表現は加算を保つので、i_M (= i·R mod p) は i_M += R_M (1 add)
// で次に進められる。M.set(i) を毎回呼ぶ必要なし。各 iter は Plantard mul 1 回
// (~3-4 cyc) + 1 add のみ。元の % p 漸化式 (~25 cyc/iter) と比べて理論上数倍速。
//
// Phase 1: fact_M[i] = mul(fact_M[i-1], i_M);  i_M += R_M
// Phase 2: cur_raw = pow(M.get(fact_M[N]), p-2); cur_M = M.set(cur_raw)  (1 modpow)
// Phase 3: inv_M[i] = mul(cur_M, fact_M[i-1]); cur_M = mul(cur_M, i_M); i_M -= R_M
// =============================================================================
struct ModInv {
 using i32 = int32_t;
 static constexpr u32 MAGIC1 = 1000000;
 static constexpr u32 MAGIC2 = 1300000;

 struct MP {
  u32 mod;
  u32 r2;
  u64 iv;
  constexpr MP(u32 m): mod(m), r2(u32(-u128(m) % m)), iv(inv_(m)) {}
  static constexpr u64 inv_(u64 n, int e = 6, u64 x = 1) {
   return e ? inv_(n, e - 1, x * (2 - x * n)) : x;
  }
  constexpr u32 reduce(u64 w) const { return u32((u128((w * iv) | u32(-1)) * mod) >> 64); }
  constexpr u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
  constexpr u32 set(u32 n) const { return mul(n, r2); }
  constexpr u32 get(u32 n) const { return reduce(n); }
 };

 static u32 pow_mod(u64 a, u64 e, u32 p) {
  u64 r = 1;
  while (e) { if (e & 1) r = r * a % p; a = a * a % p; e >>= 1; }
  return u32(r);
 }

 static vector<u32> run(u32 p, const vector<u32>& qs) {
  MP M(p);
  const u32 R_M = M.set(1);  // Montgomery 形の 1

  // Phase 1: fact_M[i] = i! · R mod p
  vector<u32> inv_M(MAGIC2 + 1);
  inv_M[0] = R_M;  // 0! = 1
  u32 i_M = R_M;
  for (u32 i = 1; i <= MAGIC2; ++i) {
   inv_M[i] = M.mul(inv_M[i - 1], i_M);
   // i_M = (i+1) · R; M.reduce が [0, mod] を返す可能性あり、p で wrap
   i_M += R_M; if (i_M >= p) i_M -= p;
  }

  // Phase 2: cur_M = ((MAGIC2)!)^{-1} · R
  u32 fact_N_raw = M.get(inv_M[MAGIC2]);
  if (fact_N_raw >= p) fact_N_raw -= p;
  u32 cur_M = M.set(pow_mod(fact_N_raw, u64(p) - 2, p));

  // Phase 3: i 降順で inv_M[i] = cur_M · (i-1)!_M
  // i_M は Phase 1 終了時に (MAGIC2+1)·R。これを R 引きつつ MAGIC2..1 を回す。
  i_M -= R_M; if (i_M >= p) i_M += p;  // 念のため。素直には i_M -= R_M で (MAGIC2)·R
  for (u32 i = MAGIC2; i >= 1; --i) {
   u32 inv_i_M = M.mul(cur_M, inv_M[i - 1]);
   cur_M = M.mul(cur_M, i_M);
   inv_M[i] = inv_i_M;
   // i_M -= R_M (= (i-1)·R)
   if (i_M < R_M) i_M += p - R_M; else i_M -= R_M;
  }
  inv_M[0] = 0;

  // ---- Farey table 構築 ----
  vector<u32> farey_lookup(MAGIC1, 0);
  auto farey_rec = [&](auto& self, u32 f1, u32 f2, u32 x, u32 y) -> void {
   u32 f3 = f1 + f2;
   u32 lo = (((u64) p * (f3 >> 16) - MAGIC2) * MAGIC1 - 1) / ((u64) p * (f3 & 0xffff)) + 1;
   u32 hi = (((u64) p * (f3 >> 16) + MAGIC2) * MAGIC1) / ((u64) p * (f3 & 0xffff));
   lo = std::max(lo, x); hi = std::min(hi, y);
   if (x < lo) self(self, f1, f3, x, lo);
   std::fill(farey_lookup.begin() + lo, farey_lookup.begin() + hi, f3);
   if (hi < y) self(self, f3, f2, hi, y);
  };
  const u32 first_x = u64(MAGIC2) * MAGIC1 / p;
  const u32 first_y = (u64(p - MAGIC2) * MAGIC1 - 1) / (p * 2) + 1;
  std::fill(farey_lookup.begin(), farey_lookup.begin() + first_x, 1u);
  farey_rec(farey_rec, 1u, 0x10002u, first_x, first_y);
  std::fill(farey_lookup.begin() + first_y, farey_lookup.begin() + MAGIC1/2, 0x10002u);
  for (u32 i = MAGIC1/2; i < MAGIC1; ++i) {
   farey_lookup[i] = (farey_lookup[MAGIC1 - 1 - i] * 0xffff0001u ^ 0xffff0000u) + 0x10000u;
  }

  // ---- per-query ----
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 x : qs) {
   if (x == 0) { ans.push_back(u32(-1)); continue; }
   if (x == 1) { ans.push_back(1); continue; }
   const u32 bucket = u64(x) * MAGIC1 / p;
   const u32 frac = farey_lookup[bucket];
   const u32 den = frac & 0xffff;
   const i32 num = (i32) (den * x - (frac >> 16) * p);
   const u32 anum = (u32) std::abs(num);
   // lazy reduction: 符号反転 (p-r) は [0, mod] → [0, mod] を保つので
   // 最終正規化 1 回で済ます
   u32 r = M.reduce(u64(den) * inv_M[anum]);
   if (num < 0) r = p - r;
   if (r >= p) r -= p;
   ans.push_back(r);
  }
  return ans;
 }
};
