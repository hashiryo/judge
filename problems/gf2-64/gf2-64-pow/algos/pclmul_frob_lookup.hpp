#pragma once
// pclmul_frob_basis.hpp の発展形 (sparse e 特化版):
//
//   α^e = ⊕(× over set bits k of e) α^{2^k}
//
// frob_basis は α^{2^k} を全 k=0..63 chain で先に生成する (sq×64) が、
// 本実装は **set bit が立っている k のみ** を byte table 1 lookup で直接生成する。
// → popcount(e) が小さい sparse e で大きな高速化。
//
// メモリ: FROB_BYTE[64][8][256] × 8 byte = 1 MiB (L2 fit)。
// 各 k で 8 lookup + 7 XOR ≈ 10 cycle、popcount mul と合わせ per query
//   popcount(e) × (10 + 5) cycle。
// frob_basis (固定 64 sq + popcount mul) と比べ:
//   popcount=1   本:    15 cycle / frob_basis: 320 + 5 = 325 cycle (~20× 速)
//   popcount=8   本:   120 cycle / frob_basis: 320 + 40 = 360 cycle (~3× 速)
//   popcount=32  本:   480 cycle / frob_basis: 320 + 160 = 480 cycle (≈ 同等)
//   popcount=64  本:   960 cycle / frob_basis: 320 + 320 = 640 cycle (~1.5× 遅)
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
namespace gf2_64_pow_frob_lookup {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
// FROB_BYTE[k][p][b] = (b << (8 p))^{2^k}  (k=0 は identity)
inline u64 FROB_BYTE[64][8][256];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 // k = 0: identity
 for(int p= 0; p < 8; ++p)
  for(int b= 0; b < 256; ++b) FROB_BYTE[0][p][b]= u64(b) << (8 * p);
 // k > 0: sq applied to k-1
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
 // set bit k のみ a^{2^k} を生成
 u64 selected[64];
 int n= 0;
 for(int k= 0; k < 64; ++k) {
  if((e >> k) & 1) selected[n++]= apply_frob(k, a);
 }
 // Tree mul (depth log2(popcount))
 while(n > 1) {
  int new_n= 0, i= 0;
  for(; i + 1 < n; i+= 2) selected[new_n++]= mul(selected[i], selected[i + 1]);
  if(i < n) selected[new_n++]= selected[i];
  n= new_n;
 }
 return selected[0];
}
}  // namespace gf2_64_pow_frob_lookup
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_frob_lookup::init_tables;
  using gf2_64_pow_frob_lookup::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
