#pragma once
// 指数 e を 16-bit ごとに 4 chunk に分けて並列処理する pow バリアント。
//
//   α^e = α^{e_0} · (α^{2^{16}})^{e_1} · (α^{2^{32}})^{e_2} · (α^{2^{48}})^{e_3}
//          where e = e_0 + e_1·2^{16} + e_2·2^{32} + e_3·2^{48}, e_i ∈ [0, 65535]
//
// 4 つの sub-pow は全て独立 (異なる base, 異なる exponent) → OoO で並列実行され
// 単一 pow の critical path が短くなる。各 sub-pow は 16-bit binary exp
// (16 sqs + ~8 muls)、frob16 byte table で base を 3 個追加生成。
//
// 期待: 単一 pow の latency は ~1/3 に短縮、ただし throughput は同等程度。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_split16 {
using gf2_64_pclmul::frob16, gf2_64_pclmul::frob32, gf2_64_pclmul::frob48;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
// 16-bit exponent の binary exp
u64 pow16(u64 a, u16 e) {
 u64 result= 1;
 while(e) {
  if(e & 1) result= mul(result, a);
  a= sq(a);
  e>>= 1;
 }
 return result;
}
u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;

 // 4 並列 sub-pow (完全に独立、OoO が並列化)
 const u64 r0= pow16(a, e);
 const u64 r1= pow16(frob16(a), e >> 16);
 const u64 r2= pow16(frob32(a), e >> 32);
 const u64 r3= pow16(frob48(a), e >> 48);
 // Tree mul (depth log4 = 2)
 return mul(mul(r0, r1), mul(r2, r3));
}
}  // namespace gf2_64_pow_split16
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_split16::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
