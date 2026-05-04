#pragma once
// sqrt(a) = a^{2^63} は GF(2)-線型写像 (Frobenius の冪)。
// 64×64 GF(2) 行列 → 8 byte tables で 1 inv あたり 8 byte lookups + XOR。
//
// 既存の pclmul.hpp は pow(a, 2^63) で 64 回二乗を回す = 64 PCLMUL+reduce ≈ 200 cycles。
// 本バリアントは ~10 cycles で完結 → 20× 速い見込み。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "../../_shared/_common.hpp"
#include "../../_shared/pclmul_core.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#include <immintrin.h>
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif
namespace gf2_64_sqrt_linear {
inline u64 SQRT_BYTE[8][256];
inline bool inited= false;
[[gnu::target("pclmul")]] void init_table() {
 if(inited) return;
 inited= true;
 // 列 j の値: (x^j)^{2^63} (= sqrt of x^j).
 // sqrt = sq^63 (= 64 sqs だと a^{2^64} = a に戻るので 63 回が正解)
 u64 col[64];
 for(int j= 0; j < 64; ++j) {
  u64 v= u64(1) << j;
  for(int k= 0; k < 63; ++k) v= gf2_64_pclmul::sq(v);
  col[j]= v;
 }
 for(int p= 0; p < 8; ++p) {
  for(int b= 0; b < 256; ++b) {
   u64 v= 0;
   for(int bit= 0; bit < 8; ++bit) {
    if((b >> bit) & 1) v^= col[p * 8 + bit];
   }
   SQRT_BYTE[p][b]= v;
  }
 }
}
[[gnu::always_inline]] inline u64 sqrt_lin(u64 a) { return SQRT_BYTE[0][u8(a)] ^ SQRT_BYTE[1][u8(a >> 8)] ^ SQRT_BYTE[2][u8(a >> 16)] ^ SQRT_BYTE[3][u8(a >> 24)] ^ SQRT_BYTE[4][u8(a >> 32)] ^ SQRT_BYTE[5][u8(a >> 40)] ^ SQRT_BYTE[6][u8(a >> 48)] ^ SQRT_BYTE[7][u8(a >> 56)]; }
}  // namespace gf2_64_sqrt_linear
struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sqrt_linear::init_table;
  using gf2_64_sqrt_linear::sqrt_lin;
  init_table();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= sqrt_lin(as[i]);
  return ans;
 }
};
