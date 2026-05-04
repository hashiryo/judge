#pragma once
// pclmul_frob_lookup.hpp の変種: __builtin_ctzll ループ + 逐次 mul (tree なし)。
// _ctz と _serial の両方の変更を取り込んだ最小コード版。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_frob_lookup_ctz_serial {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
inline u64 FROB_BYTE[64][8][256];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 for(int p= 0; p < 8; ++p)
  for(int b= 0; b < 256; ++b) FROB_BYTE[0][p][b]= u64(b) << (8 * p);
 for(int k= 1; k < 64; ++k)
  for(int p= 0; p < 8; ++p)
   for(int b= 0; b < 256; ++b) FROB_BYTE[k][p][b]= sq(FROB_BYTE[k - 1][p][b]);
}
inline u64 apply_frob(int k, u64 a) {
 const auto& t= FROB_BYTE[k];
 return t[0][u8(a)] ^ t[1][u8(a >> 8)] ^ t[2][u8(a >> 16)] ^ t[3][u8(a >> 24)] ^ t[4][u8(a >> 32)] ^ t[5][u8(a >> 40)] ^ t[6][u8(a >> 48)] ^ t[7][u8(a >> 56)];
}
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 u64 m= e;
 int k0= __builtin_ctzll(m);
 m&= m - 1;
 u64 acc= apply_frob(k0, a);
 while(m) {
  int k= __builtin_ctzll(m);
  m&= m - 1;
  acc= mul(acc, apply_frob(k, a));
 }
 return acc;
}
}  // namespace gf2_64_pow_frob_lookup_ctz_serial
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_frob_lookup_ctz_serial::init_tables;
  using gf2_64_pow_frob_lookup_ctz_serial::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
