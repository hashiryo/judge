#pragma once
#include "_common.hpp"
// 標準的な Baby-step Giant-step。
// 全 query で同じ baby step テーブル (g^j for j=0..m-1) を共有。
// per query: m ≈ √p ≈ 31623 iterations。p ≈ 10^9 query 数 n に対し O(n √p)。
struct DLog {
 static u32 mod_pow(u32 b, u64 e, u32 mod) {
  u64 r = 1, x = b;
  while (e) { if (e & 1) r = r * x % mod; x = x * x % mod; e >>= 1; }
  return u32(r);
 }
 static vector<u32> run(u32 p, u32 g, const vector<u32>& qs) {
  // m = ceil(sqrt(p-1))
  u32 m = 1;
  while (u64(m) * m < p - 1) ++m;
  // 開放アドレス hash table (load 0.25)
  u64 cap = 8;
  while (cap < 4ull * m) cap *= 2;
  vector<u32> keys(cap, 0);
  vector<u32> values(cap, u32(-1));
  u64 mask = cap - 1;
  u64 cur = 1;
  for (u32 j = 0; j < m; ++j) {
   u64 h = (u64(u32(cur)) * 0x9E3779B9u) & mask;
   while (values[h] != u32(-1)) h = (h + 1) & mask;
   keys[h] = u32(cur);
   values[h] = j;
   cur = cur * g % p;
  }
  const u32 inv_step = mod_pow(g, p - 1 - m, p);  // = g^{-m}
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) {
   if (a == 0) { ans.push_back(u32(-1)); continue; }
   u64 t = a;
   bool found = false;
   for (u32 i = 0; i < m && !found; ++i) {
    u64 h = (u64(u32(t)) * 0x9E3779B9u) & mask;
    while (values[h] != u32(-1)) {
     if (keys[h] == t) {
      u64 res = u64(i) * m + values[h];
      if (res < p - 1) { ans.push_back(u32(res)); found = true; break; }
     }
     h = (h + 1) & mask;
    }
    if (!found) t = t * inv_step % p;
   }
   if (!found) ans.push_back(u32(-1));
  }
  return ans;
 }
};
