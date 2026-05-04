#pragma once
// Itoh-Tsujii 逆元。
//
// 素朴な Fermat 逆元 a^{2^64 - 2} は 63 muls + 64 sqs (~127 ops) かかるが、
// 中間値 T_k = a^{2^k - 1} を再利用する addition chain で 10 muls + 63 sqs (~73 ops)
// に削減できる。 理論値 1.7× 高速化。
//
// mul / sq は _shared/pclmul_core.hpp の current best を使用 (sq は PDEP-based)。
//
// チェーン構造:
//   main: T_1, T_2, T_4, T_8, T_16, T_32 を doubling で構築 (T_{2k} = T_k · T_k^{2^k})
//         5 muls, 31 sqs
//   side: T_63 を mixed-radix 展開 (63 = 32+16+8+4+2+1) で構築
//         acc = T_32 → ((acc^{2^16})·T_16)^{2^8}·T_8 → ... → acc·T_1
//         5 muls, 31 sqs
//   final: a^{2^64 - 2} = T_63^2 (1 sq)
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pclmul_itoh_tsujii {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob8;
using gf2_64_pclmul::frob16;
u64 inv(u64 a) {
 // main chain: T_k = a^{2^k - 1} for k = 1, 2, 4, 8, 16, 32
 // frobK(_) = sq×K の byte-table 版 (_shared/frob.hpp)
 const u64 T1= a;
 const u64 T2= mul(T1, sq(T1));         //  1 sq + 1 mul
 const u64 T4= mul(T2, frob2(T2));      //  frob2 + 1 mul
 const u64 T8= mul(T4, frob4(T4));      //  frob4 + 1 mul
 const u64 T16= mul(T8, frob8(T8));     //  frob8 + 1 mul
 const u64 T32= mul(T16, frob16(T16));  // frob16 + 1 mul
 // side chain: T_63 = a^{2^63 - 1} を 63 = 32+16+8+4+2+1 で組み立てる
 u64 acc= mul(frob16(T32), T16);  // a^{2^48 - 1}
 acc= mul(frob8(acc), T8);        // a^{2^56 - 1}
 acc= mul(frob4(acc), T4);        // a^{2^60 - 1}
 acc= mul(frob2(acc), T2);        // a^{2^62 - 1}
 acc= mul(sq(acc), T1);           // a^{2^63 - 1}
 return sq(acc);                  // a^{2^64 - 2} = a^{-1}
}
}  // namespace gf2_64_pclmul_itoh_tsujii
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_itoh_tsujii::inv;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
