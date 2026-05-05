#pragma once
// PCLMUL + 4-bit window + Frobenius byte table。
//
// pow(a, e) は a と e 両方が変数なので linear_map で完全 byte 化はできない。
// しかし「multi-step squaring (Frobenius)」部分は GF(2)-線型なので byte table
// 化できる。これにより:
//   - 64 個の sq をシーケンシャルに繋ぐ部分が消える (64 × 2 cycle = 128 cycle)
//   - 代わりに 15 個の frob_4 (= a → a^{16}) を byte table で適用
//     (各 8 lookups, ~10 cycles, 15 回の chain で 150 cycles)
//
// 加えて 4-bit window で precompute (a^0..a^15) しておけば、e の各 nibble を
// 1 mul で消化できる。
//
// 期待: pclmul_pdep より速い、特に arm-clang (PDEP fallback が遅い) で大きく勝つ。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_byte_window {
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
u64 pow(u64 a, u64 e) {
 if(e == 0) return 1;
 // Precompute T[i] = a^i for i = 0..15 (14 muls)
 u64 T[16]= {1, a};
 for(int i= 2; i < 16; ++i) T[i]= mul(T[i - 1], a);

 // Find top nonzero nibble of e
 int top= 15 - (__builtin_clzll(e) >> 2);

 u64 acc= T[(e >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob4(acc);
  u16 chunk= (e >> (4 * i)) & 0xF;
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
}  // namespace gf2_64_pow_byte_window
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_byte_window::pow;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
