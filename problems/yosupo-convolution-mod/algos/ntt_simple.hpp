#pragma once
#include "_common.hpp"
// 教科書的反復 NTT。Montgomery reduction なしの「u64 で受けて u32 mod する」型。
// 998244353 = 119 * 2^23 + 1, primitive root g = 3。
struct Conv {
 static u32 pow_mod(u32 a, u64 e) {
  u64 r= 1, x= a;
  while (e) {
   if (e & 1) r= r * x % MOD;
   x= x * x % MOD;
   e>>= 1;
  }
  return (u32)r;
 }
 static void ntt(vector<u32>& a, bool inv) {
  int n= (int)a.size();
  for (int i= 1, j= 0; i < n; ++i) {
   int bit= n >> 1;
   for (; j & bit; bit>>= 1) j^= bit;
   j^= bit;
   if (i < j) swap(a[i], a[j]);
  }
  for (int len= 2; len <= n; len<<= 1) {
   u32 w= inv ? pow_mod(pow_mod(3, MOD - 1 - (MOD - 1) / len), 1) : pow_mod(3, (MOD - 1) / len);
   for (int i= 0; i < n; i+= len) {
    u64 wn= 1;
    for (int j= 0; j < len / 2; ++j) {
     u32 u= a[i + j], v= (u32)(wn * a[i + j + len / 2] % MOD);
     a[i + j]= u + v < MOD ? u + v : u + v - MOD;
     a[i + j + len / 2]= u >= v ? u - v : u + MOD - v;
     wn= wn * w % MOD;
    }
   }
  }
  if (inv) {
   u64 ninv= pow_mod(n, MOD - 2);
   for (auto& x: a) x= (u32)(x * ninv % MOD);
  }
 }
 static vector<u32> run(const vector<u32>& a, const vector<u32>& b) {
  int n= (int)a.size(), m= (int)b.size();
  if (!n || !m) return {};
  int sz= n + m - 1, len= 1;
  while (len < sz) len<<= 1;
  vector<u32> fa(a.begin(), a.end()), fb(b.begin(), b.end());
  fa.resize(len, 0);
  fb.resize(len, 0);
  ntt(fa, false);
  ntt(fb, false);
  for (int i= 0; i < len; ++i) fa[i]= (u32)((u64)fa[i] * fb[i] % MOD);
  ntt(fa, true);
  fa.resize(sz);
  return fa;
 }
};
