#pragma once
#include "_common.hpp"
// Plantard 32-bit + base 側を preconditioned 形 (base * iv mod 2^64) で持って
// pow ループ内の mul で 1 mul/call 削減する変種。
//
// 標準 Plantard mul(l, r) は:
//   1. l * r → u64 (1 mul)
//   2. * iv → low u64 (1 mul)
//   3. | u32(-1) → u128
//   4. * mod → u128 (~2 mul)
//   合計 4 mul/call。
//
// preconditioned 形 r_pre = r * iv mod 2^64 を保持すると:
//   mul_pre(l, r_pre):
//   1. l * r_pre → low u64 (1 mul) ← step 1+2 を 1 mul に圧縮
//   2. | u32(-1)
//   3. * mod (~2 mul)
//   合計 3 mul/call。1 mul 削減。
//
// pow ループ全体では: r 更新 (mul_pre) + base 更新 (mul_pre + to_pre) で
// 標準 8 mul/iter → 7 mul/iter 程度。
//
// 参考: judge/problems/modpow-test/runtime/runtime-62/algos/plantard64_pre.hpp
struct ModInv {
 struct MP {  // mod < 2^32/phi
  u32 mod;
  u32 r2;
  u64 iv;
  constexpr MP(u32 m): mod(m), r2(u32(-u128(m) % m)), iv(inv_(m)) {}
  static constexpr u64 inv_(u64 n, int e = 6, u64 x = 1) {
   return e ? inv_(n, e - 1, x * (2 - x * n)) : x;
  }
  // 標準 reduce
  constexpr u32 reduce(u64 w) const { return u32((u128((w * iv) | u32(-1)) * mod) >> 64); }
  constexpr u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
  constexpr u32 set(u32 n) const { return mul(n, r2); }
  constexpr u32 get(u32 n) const { return reduce(n); }
  // preconditioned: r_pre = r * iv mod 2^64
  // mul_pre(l, r_pre): l * r_pre は (l * r) * iv mod 2^64 と等価
  inline u64 to_pre(u32 r) const { return u64(r) * iv; }
  inline u32 mul_pre(u32 l, u64 r_pre) const {
   const u64 t = u64(l) * r_pre;  // = (l * r * iv) mod 2^64
   return u32((u128(t | u32(-1)) * mod) >> 64);
  }
  // Fermat 用 pow: r 更新は mul_pre 経由、base 更新は base × base + to_pre 再計算
  inline u32 pow(u32 base, u32 e) const {
   u32 r = set(1);
   if (!e) return r;
   u64 base_pre = to_pre(base);
   for (;;) {
    if (e & 1) r = mul_pre(r, base_pre);
    if (!(e >>= 1)) break;
    base = mul_pre(base, base_pre);  // base ← base * base via preconditioned
    base_pre = to_pre(base);          // 新 base 用の preconditioned 値
   }
   return r;
  }
 };
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  MP M(p);
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) {
   if (a == 0) { ans.push_back(u32(-1)); continue; }
   const u32 a_m = M.set(a);
   const u32 inv_m = M.pow(a_m, p - 2);
   ans.push_back(M.get(inv_m));
  }
  return ans;
 }
};
