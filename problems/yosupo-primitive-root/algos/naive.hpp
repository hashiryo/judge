#pragma once
#include "_common.hpp"
// 素朴: g = 2, 3, ... を順に試して order(g) == p-1 になる最初の g を返す。
// p-1 の因数分解は試し割り (~ 1e9 まで)。large p では TLE 必至だが、
// 比較対象として置いておく。

namespace primitive_root_naive {
inline u64 mulmod(u64 a, u64 b, u64 m) { return (__uint128_t)a * b % m; }
inline u64 powmod(u64 a, u64 e, u64 m) {
 u64 r = 1 % m;
 a %= m;
 while (e) {
  if (e & 1) r = mulmod(r, a, m);
  a = mulmod(a, a, m);
  e >>= 1;
 }
 return r;
}
inline vector<u64> factor_pminus1_naive(u64 n) {
 // 素朴 trial division。p-1 の素因数を unique で返す。
 vector<u64> ps;
 for (u64 p = 2; p * p <= n && p < 1000000; ++p) {
  if (n % p == 0) {
   ps.push_back(p);
   while (n % p == 0) n /= p;
  }
 }
 if (n > 1) ps.push_back(n);
 return ps;
}
inline u64 primitive_root(u64 p) {
 if (p == 2) return 1;
 auto ps = factor_pminus1_naive(p - 1);
 for (u64 g = 2;; ++g) {
  bool ok = true;
  for (u64 q : ps) {
   if (powmod(g, (p - 1) / q, p) == 1) { ok = false; break; }
  }
  if (ok) return g;
 }
}
} // namespace

struct PrimitiveRoot {
 static vector<u64> run(const vector<u64>& qs) {
  vector<u64> ans;
  ans.reserve(qs.size());
  for (auto p : qs) ans.push_back(primitive_root_naive::primitive_root(p));
  return ans;
 }
};
