#pragma once
// Norm-based inv に「Frobenius byte table を 3 段 (frob16/32/48) 持つ」並列化。
//
// pclmul_norm.hpp は frob16 byte table を 3 回チェイン適用するため依存チェーン上
// 24 lookups が直列。frob32 / frob48 の byte table も持つことで β1, β2, β3 を
// 独立に計算でき、CPU の OoO で並列実行されることを期待。
//
// メモリ: 3 × 16 KiB = 48 KiB byte tables (L1 ぎりぎり)。
// 計算量: 同じく 24 lookups 合計、ただし依存チェーン解消。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"

namespace gf2_64_pclmul_norm_parallel {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob48;

// F_{2^16} 内の inv: N^{-1} = N^{2^16 - 2}, addition chain on 15 = 8+4+2+1
u64 inv_in_f16(u64 N) {
 const u64 T1= N;
 const u64 T2= mul(T1, sq(T1));
 const u64 T4= mul(T2, frob2(T2));
 const u64 T8= mul(T4, frob4(T4));
 u64 acc= mul(frob4(T8), T4);
 acc= mul(frob2(acc), T2);
 acc= mul(sq(acc), T1);
 return sq(acc);
}
u64 inv(u64 a) {
 // β1, β2, β3 を独立 byte table で並列計算 (依存なし → ILP 効くはず)
 const u64 b1= frob16(a);
 const u64 b2= frob32(a);
 const u64 b3= frob48(a);
 const u64 beta= mul(mul(b1, b2), b3);
 const u64 N= mul(a, beta);
 return mul(beta, inv_in_f16(N));
}
}  // namespace gf2_64_pclmul_norm_parallel
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_parallel::inv;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
