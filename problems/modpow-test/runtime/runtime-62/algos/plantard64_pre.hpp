#pragma once
#include "_common.hpp"
// Plantard 乗算で base 側を preconditioned 形 (base * iv mod 2^128) で持って
// 乗算 1 回あたり 2 mul 削る変種。標準の plantard64.hpp と数値的には等価。
//
// 標準 mul(l, r) = reduce(u128(l) * r) は内部で
//   t = (l*r * iv) mod 2^128                 (u128*u128 low128 = 3 muls)
//   q = t >> 64
//   return (u128(q+1) * mod) >> 64            (1 mul)
// で 4 muls + 引数準備の 1 mul = 5 muls (u64*u64) 程度。
//
// base_pre = base * iv mod 2^128 を持つと
//   t = (l * base_pre) mod 2^128              (u64*u128 low128 = 2 muls)
//   q, return                                  (1 mul)
// で 3 muls。pow ループの「r *= base」側でこれを使う。
// squaring 側は base_new = mul(base, base) を出した後に
// base_new_pre = base_new * iv mod 2^128 (2 muls) を更新するので、合計は
// 同じ 5 muls だが、r *= base 部分で 2 mul / call 削れる。
struct MP {  // mod < 2^64/phi
 u64 mod;
 constexpr MP(): mod(0), r2(0), iv(0) {}
 constexpr MP(u64 m): mod(m), r2(R2(m)), iv(inv(m)) {}
 constexpr inline u64 mul(u64 l, u64 r) const { return reduce(u128(l) * r); }
 constexpr inline u64 set(u64 n) const { return mul(n, r2); }
 constexpr inline u64 get(u64 n) const { return n= reduce(n), n >= mod ? n - mod : n; }
 constexpr inline u64 norm(u64 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u64 plus(u64 l, u64 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u64 diff(u64 l, u64 r) const { return l-= r, l >> 63 ? l + mod : l; }
 inline u64 pow(u64 base, u64 e) const {
  u64 r = set(1);
  if (!e) return r;
  u128 base_pre = to_pre(base);
  for (;;) {
   if (e & 1) r = mul_pre(r, base_pre);
   if (!(e >>= 1)) break;
   base = mul_pre(base, base_pre);
   base_pre = to_pre(base);
  }
  return r;
 }
private:
 u64 r2;
 u128 iv;
 // l * r_pre mod 2^128 → reduce した値 (∈ [0, mod])。
 // l は u64 なので u128(l)*r_pre は実質 2 mul (hi(l)*... が 0 で潰れる)。
 inline u64 mul_pre(u64 l,const u128& r_pre) const {
  return (u128(u64((l * r_pre) >> 64) + 1) * mod) >> 64;
 }
 // x * iv mod 2^128。x ∈ [0, mod) なら preconditioned 形。
 inline u128 to_pre(u64 x) const { return x * iv; }
 static constexpr u128 inv(u128 n, int e= 7, u128 x= 1) { return e ? inv(n, e - 1, x * (2 - x * n)) : x; }
 constexpr inline u64 reduce(const u128& w) const { return (u128(u64((w * iv) >> 64) + 1) * mod) >> 64; }
 static constexpr u64 R2(u64 n) {
  u64 r= -u128(n) % n;
  return u128(r) * r % n;
 }
};
