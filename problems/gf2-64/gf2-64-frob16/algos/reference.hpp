#pragma once
#include "../../_shared/_common.hpp"
// 素朴 reference: bit-by-bit sq を 16 回繰り返し。
struct GF2_64Op {
 static u64 sq_naive(u64 a) {
  u64 lo= 0, hi= 0;
  for(int i= 0; i < 64; ++i) {
   if((a >> i) & 1) {
    int j= 2 * i;
    if(j < 64) lo^= u64(1) << j;
    else hi^= u64(1) << (j - 64);
   }
  }
  for(int i= 63; i >= 0; --i) {
   if((hi >> i) & 1) {
    hi^= u64(1) << i;
    lo^= IRRED_LOW << i;
    if(i > 0) hi^= IRRED_LOW >> (64 - i);
   }
  }
  return lo;
 }
 static u64 frob16(u64 a) {
  for(int i= 0; i < 16; ++i) a= sq_naive(a);
  return a;
 }
 static vector<u64> run(const vector<u64>& as) {
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= frob16(as[i]);
  return ans;
 }
};
