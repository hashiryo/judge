#pragma once
// Norm-based inversion (subfield reduction) on poly basis.
//
// 塔 F_{2^64} ⊃ F_{2^16} を活用。norm を使うと F_{2^64}-inv を F_{2^16}-inv に reduce できる:
//
//   N(α) = α · α^{2^16} · α^{2^32} · α^{2^48}  ∈ F_{2^16}  (Galois 共役の積は subfield)
//   α^{-1} = (α^{2^16} · α^{2^32} · α^{2^48}) · N(α)^{-1}
//
// 各 Frobenius α^{2^16} は GF(2)-線型なので 64×64 行列 → 8 byte tables (16 KiB) で
// 8 回 byte lookup に置換可能。16 PCLMUL-square より少ない命令数で済む。
//
// N は F_{2^16} subfield に住む。inv は F_{2^64} 上の Itoh-Tsujii で計算するが、
// N^{2^16 - 1} = 1 なので a^{-1} = a^{2^16 - 2} で済み、F_{2^15 - 1} までの addition chain
// で十分 (= 6 muls + 15 sqs)。
//
// 期待 ops: 10 muls + 15 sqs + 24 byte lookups
// (cf. pclmul_itoh_tsujii: 10 muls + 63 sqs)
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"

namespace gf2_64_pclmul_norm {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob16;

// F_{2^16} 内の inv: N^{-1} = N^{2^16 - 2} = (N^{2^15 - 1})^2
// 15 = 8+4+2+1 で addition chain (Itoh-Tsujii 縮小版)
u64 inv_in_f16(u64 N) {
 const u64 T1= N;
 const u64 T2= mul(T1, sq(T1));     // a^{2^2 - 1}
 const u64 T4= mul(T2, frob2(T2));  // a^{2^4 - 1}
 const u64 T8= mul(T4, frob4(T4));  // a^{2^8 - 1}
 u64 acc= mul(frob4(T8), T4);       // a^{2^12 - 1}
 acc= mul(frob2(acc), T2);          // a^{2^14 - 1}
 acc= mul(sq(acc), T1);             // a^{2^15 - 1}
 return sq(acc);                    // a^{2^16 - 2} = a^{-1}
}
u64 inv(u64 a) {
 // β1, β2, β3 = α^{2^16}, α^{2^32}, α^{2^48}, frob16 で順に
 const u64 b1= frob16(a);
 const u64 b2= frob16(b1);
 const u64 b3= frob16(b2);
 const u64 beta= mul(mul(b1, b2), b3);
 const u64 N= mul(a, beta);
 return mul(beta, inv_in_f16(N));
}
}  // namespace gf2_64_pclmul_norm
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm::inv;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
