#pragma once
// PCLMUL + Itoh-Tsujii + PDEP 二乗。
//
// GF(2) の二乗は (Σ a_i x^i)^2 = Σ a_i x^{2i} と cross term が消える (Frobenius)。
// したがって 64-bit 入力 → 128-bit 出力は単なる「bit を 1 つおきに spread する」
// 操作で、PDEP (BMI2) 1 命令で書ける。x86-64-v3 (= BMI2 あり) なら PCLMUL より
// 安く、Itoh-Tsujii の 63 sqs を大幅に短縮できる。
//
// ARM では PDEP が無いので bit-twiddling fallback (5 SHIFT + 5 AND)。
// それでも PCLMUL + reduce より軽い。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>  // _pdep_u64
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif
namespace gf2_64_pclmul_itoh_tsujii_pdep {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
[[gnu::target("pclmul")]] u64 frob(u64 x, int k) {
 for(int i= 0; i < k; ++i) x= sq(x);
 return x;
}
[[gnu::target("pclmul")]] u64 inv(u64 a) {
 // main chain: T_k = a^{2^k - 1} for k = 1, 2, 4, 8, 16, 32
 const u64 T1= a;
 const u64 T2= mul(T1, frob(T1, 1));
 const u64 T4= mul(T2, frob(T2, 2));
 const u64 T8= mul(T4, frob(T4, 4));
 const u64 T16= mul(T8, frob(T8, 8));
 const u64 T32= mul(T16, frob(T16, 16));
 // side chain: T_63 を 63 = 32+16+8+4+2+1 で組み立てる
 u64 acc= mul(frob(T32, 16), T16);  // a^{2^48 - 1}
 acc= mul(frob(acc, 8), T8);        // a^{2^56 - 1}
 acc= mul(frob(acc, 4), T4);        // a^{2^60 - 1}
 acc= mul(frob(acc, 2), T2);        // a^{2^62 - 1}
 acc= mul(frob(acc, 1), T1);        // a^{2^63 - 1}
 return sq(acc);                    // a^{2^64 - 2} = a^{-1}
}
}  // namespace gf2_64_pclmul_itoh_tsujii_pdep
struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_itoh_tsujii_pdep::inv;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
