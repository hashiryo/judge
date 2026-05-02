#pragma once
// PCLMUL + Itoh-Tsujii 逆元。
//
// 素朴な Fermat 逆元 a^{2^64 - 2} は 63 muls + 64 sqs (~127 ops) かかるが、
// 中間値 T_k = a^{2^k - 1} を再利用する addition chain で 10 muls + 63 sqs (~73 ops)
// に削減できる。理論値 1.7× 高速化。
//
// チェーン構造:
//   main: T_1, T_2, T_4, T_8, T_16, T_32 を doubling で構築 (T_{2k} = T_k · T_k^{2^k})
//         5 muls, 31 sqs
//   side: T_63 を mixed-radix 展開 (63 = 32+16+8+4+2+1) で構築
//         acc = T_32 → ((acc^{2^16})·T_16)^{2^8}·T_8 → ... → acc·T_1
//         5 muls, 31 sqs
//   final: a^{2^64 - 2} = T_63^2 (1 sq)
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#define PCLMUL_RUN [[gnu::target("pclmul")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_pclmul_itoh_tsujii {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// k 回 Frobenius (= 二乗)。k=0 は何もしない。
PCLMUL_FN u64 frob(u64 x, int k) {
 for (int i = 0; i < k; ++i) x = sq(x);
 return x;
}

PCLMUL_FN u64 inv(u64 a) {
 // main chain: T_k = a^{2^k - 1} for k = 1, 2, 4, 8, 16, 32
 const u64 T1 = a;
 const u64 T2 = mul(T1, frob(T1, 1));   //  1 sq + 1 mul
 const u64 T4 = mul(T2, frob(T2, 2));   //  2 sqs + 1 mul
 const u64 T8 = mul(T4, frob(T4, 4));   //  4 sqs + 1 mul
 const u64 T16 = mul(T8, frob(T8, 8));  //  8 sqs + 1 mul
 const u64 T32 = mul(T16, frob(T16, 16)); // 16 sqs + 1 mul
 // side chain: T_63 = a^{2^63 - 1} を 63 = 32+16+8+4+2+1 で組み立てる
 //   acc が a^{2^k - 1} のとき、acc = mul(frob(acc, j), T_j) すると a^{2^{k+j} - 1}
 u64 acc = mul(frob(T32, 16), T16);  // a^{2^48 - 1}, 16 sqs + 1 mul
 acc = mul(frob(acc,  8), T8);       // a^{2^56 - 1},  8 sqs + 1 mul
 acc = mul(frob(acc,  4), T4);       // a^{2^60 - 1},  4 sqs + 1 mul
 acc = mul(frob(acc,  2), T2);       // a^{2^62 - 1},  2 sqs + 1 mul
 acc = mul(frob(acc,  1), T1);       // a^{2^63 - 1},  1 sq  + 1 mul
 return sq(acc);                     // a^{2^64 - 2} = a^{-1}, 1 sq
}

} // namespace gf2_64_pclmul_itoh_tsujii

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_itoh_tsujii::inv;
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = mul(as[i], inv(bs[i]));
  return ans;
 }
};
