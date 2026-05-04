#pragma once
// sq(a) = mul(a, a) を PCLMUL で計算 (素朴ベースライン)。
// 自己完結 (= _shared/pclmul_core.hpp に依存しない、 sq の比較対象なので)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
namespace gf2_64_sq_pclmul_mul {
inline u64 sq(u64 a) {
 __m128i av= _mm_cvtsi64_si128(a);
 __m128i v= _mm_clmulepi64_si128(av, av, 0);
 u64 h= (u64)v[1], d= h ^ (h << 1);
 return (u64)v[0] ^ ((u8[]){0, 27, 45, 54, 90, 65, 119, 108})[h >> 60] ^ d ^ (d << 3);
}
}
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sq_pclmul_mul::sq;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sq(as[i]);
  return ans;
 }
};
