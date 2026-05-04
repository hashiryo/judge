#pragma once
// pclmul_frob_lookup.hpp の変種: tree mul を廃止して逐次 mul。
// 0..63 ループはそのまま。
//
// 直列 mul は depth = popcount(e) - 1 mul、tree は log2(popcount) mul。
// dense e (popcount ~32) では tree 5 段 vs 直列 31 段なので tree 圧勝のはずだが、
// 配列管理を省ける分のオーバーヘッド差を実機で見る目的の比較対象。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_frob_lookup_serial {
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
[[gnu::always_inline]] inline u64 apply_frob(int k, u64 a) {
 const auto& t= FROB_BYTE[k];
 return t[0][u8(a)] ^ t[1][u8(a >> 8)] ^ t[2][u8(a >> 16)] ^ t[3][u8(a >> 24)] ^ t[4][u8(a >> 32)] ^ t[5][u8(a >> 40)] ^ t[6][u8(a >> 48)] ^ t[7][u8(a >> 56)];
}
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 u64 acc= 0;
 bool started= false;
 for(int k= 0; k < 64; ++k) {
  if((e >> k) & 1) {
   const u64 v= apply_frob(k, a);
   acc= started ? mul(acc, v) : v;
   started= true;
  }
 }
 return acc;
}
}  // namespace gf2_64_pow_frob_lookup_serial
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_frob_lookup_serial::init_tables;
  using gf2_64_pow_frob_lookup_serial::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
