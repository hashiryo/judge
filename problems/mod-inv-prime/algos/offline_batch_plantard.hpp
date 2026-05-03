#pragma once
#include "_common.hpp"
// offline_batch を Plantard 32-bit で高速化。
// 全 mul (prefix product + back substitution + Fermat) を Plantard reduce 経由に。
// per query 償却 ~3 mul × 10 cycle ≈ 30 cycle (vs 標準 % p の 30 cycle × 3 = 90 cycle)
struct ModInv {
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
  constexpr u32 pow(u32 base, u32 e) const {
   u32 r = set(1);
   while (e) {
    if (e & 1) r = mul(r, base);
    base = mul(base, base);
    e >>= 1;
   }
   return r;
  }
 };
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  MP M(p);
  const int Q = (int) qs.size();
  vector<u32> ans(Q, u32(-1));
  vector<u32> prefix(Q + 1);  // Montgomery form
  prefix[0] = M.set(1);
  for (int i = 0; i < Q; ++i) {
   u32 a = qs[i];
   if (a == 0) prefix[i + 1] = prefix[i];
   else prefix[i + 1] = M.mul(prefix[i], M.set(a));
  }
  u32 total_inv = M.pow(prefix[Q], p - 2);  // Montgomery 内で Fermat
  for (int i = Q - 1; i >= 0; --i) {
   u32 a = qs[i];
   if (a == 0) { ans[i] = u32(-1); continue; }
   const u32 inv_m = M.mul(prefix[i], total_inv);
   ans[i] = M.get(inv_m);
   total_inv = M.mul(total_inv, M.set(a));
  }
  return ans;
 }
};
