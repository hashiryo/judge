#pragma once
#include "_common.hpp"
// 決定的 Miller-Rabin (u64 全域で正しい 7 base): {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}
// 試し割り (小さい素数) と Miller-Rabin の組み合わせ。
// 64-bit 乗算は __int128 で素朴に。
struct Primality {
 static u64 mulmod_(u64 a, u64 b, u64 m) {
  return (u128) a * b % m;
 }
 static u64 powmod_(u64 a, u64 d, u64 m) {
  u64 r = 1;
  while (d) { if (d & 1) r = mulmod_(r, a, m); a = mulmod_(a, a, m); d >>= 1; }
  return r;
 }
 static bool miller_rabin(u64 n, u64 a) {
  if (n % a == 0) return n == a;
  u64 d = n - 1;
  int s = __builtin_ctzll(d);
  d >>= s;
  u64 x = powmod_(a, d, n);
  if (x == 1 || x == n - 1) return true;
  for (int i = 0; i < s - 1; ++i) {
   x = mulmod_(x, x, n);
   if (x == n - 1) return true;
  }
  return false;
 }
 static bool is_prime(u64 n) {
  if (n < 2) return false;
  for (u64 p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
   if (n == p) return true;
  for (u64 p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
   if (n % p == 0) return false;
  for (u64 a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
   if (!miller_rabin(n, a)) return false;
  return true;
 }
 static vector<bool> run(const vector<u64>& qs) {
  vector<bool> ans(qs.size());
  for (size_t i = 0; i < qs.size(); ++i) ans[i] = is_prime(qs[i]);
  return ans;
 }
};
