#pragma once
#include "_common.hpp"
// 試し割り: g^k = a を k=0..p-2 で順に試す。
//   p ≈ 10^9 なので 1 query あたり O(p) で激遅。比較 baseline。
struct DLog {
 static vector<u32> run(u32 p, u32 g, const vector<u32>& qs) {
  vector<u32> ans;
  ans.reserve(qs.size());
  for (u32 a : qs) {
   if (a == 0) { ans.push_back(u32(-1)); continue; }
   if (a == 1) { ans.push_back(0); continue; }
   u64 cur = 1;
   u32 k = 0;
   while (cur != a) {
    cur = cur * g % p;
    ++k;
    if (k >= p - 1) { k = u32(-1); break; }
   }
   ans.push_back(k);
  }
  return ans;
 }
};
