#pragma once
#include "_common.hpp"
// Offline batch inverse: 全 query を先読みして O(Q + log p) で全部の逆元を求める。
//
// アイデア:
//   prefix[i] = a_1 · a_2 · ... · a_i  (modular product, O(Q) muls)
//   total_inv = prefix[Q]^{-1}          (one Fermat, O(log p))
//   for i = Q downto 1:
//     a_i_inv = prefix[i-1] · total_inv  // prefix[i-1] · (a_1·..·a_i)^{-1} = a_i^{-1}
//     total_inv *= a_i                   // 残り部分の逆元更新
//
// per-query 償却: ~3 mul + 0 div = O(1)。Fermat は 1 回だけ。
// 0 入力が混じる場合: 0 を product から除外して別管理 (= 0 扱い、戻り値 -1)。
struct ModInv {
 static u32 mod_pow(u32 b, u32 e, u32 m) {
  u64 r = 1, x = b;
  while (e) { if (e & 1) r = r * x % m; x = x * x % m; e >>= 1; }
  return u32(r);
 }
 static vector<u32> run(u32 p, const vector<u32>& qs) {
  const int Q = (int) qs.size();
  vector<u32> ans(Q, u32(-1));
  vector<u32> prefix(Q + 1);
  prefix[0] = 1;
  for (int i = 0; i < Q; ++i) {
   u32 a = qs[i];
   prefix[i + 1] = (a == 0) ? prefix[i] : u32(u64(prefix[i]) * a % p);
  }
  u32 total_inv = mod_pow(prefix[Q], p - 2, p);  // 1 度の Fermat
  for (int i = Q - 1; i >= 0; --i) {
   u32 a = qs[i];
   if (a == 0) { ans[i] = u32(-1); continue; }
   ans[i] = u32(u64(prefix[i]) * total_inv % p);
   total_inv = u32(u64(total_inv) * a % p);
  }
  return ans;
 }
};
