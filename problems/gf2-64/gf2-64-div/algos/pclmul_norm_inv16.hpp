#pragma once
// pclmul_norm_logexp.hpp の発展形:
// F_{2^16} subfield の inv を log/exp 2 lookups → 直接 inv table 1 lookup に置換。
//
// memory: 256 KiB (PW + LN) → 128 KiB (INV16 のみ) で半減 + 1 lookup 節約。
// 1 inv あたり ~12 cycle 短縮見込み。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
#include "../../_shared/basis_change.hpp"
namespace gf2_64_pclmul_norm_inv16 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
// F_{2^16} 直接 inv テーブル (nim 基底): INV16[A] = A^{-1} for A != 0, INV16[0] = 0.
inline u16 INV16[65536];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 // F_{2^16} の log/exp を一旦構築 → INV16 を作る → log/exp は捨てる
 u16 PW[65536], LN[65536];
 PW[0]= PW[65535]= 1;
 for(int i= 1; i < 65535; ++i) {
  PW[i]= u16((PW[i - 1] << 1) ^ (0x1681fu & u16(-(PW[i - 1] >= 0x8000u))));
 }
 constexpr u16 f2n[16]= {0x0001u, 0x2827u, 0x392bu, 0x8000u, 0x20fdu, 0x4d1du, 0xde4au, 0x0a17u, 0x3464u, 0xe3a9u, 0x6d8du, 0x34bcu, 0xa921u, 0xa173u, 0x0ebcu, 0x0e69u};
 for(int i= 1; i < 65535; ++i) {
  u16 x= PW[i], y= 0;
  for(; x; x&= x - 1) y^= f2n[__builtin_ctz(x)];
  PW[i]= y;
  LN[y]= u16(i);
 }
 LN[1]= 0;
 INV16[0]= 0;
 INV16[1]= 1;
 for(int v= 2; v < 65536; ++v) {
  INV16[v]= PW[(65535u - u32(LN[v])) % 65535u];
 }
}
u64 inv_in_f16_direct(u64 N_poly) {
 const u64 N_nim= gf2_64_basis::poly_to_nim(N_poly);
 const u16 inv_n16= INV16[u16(N_nim)];  // 1 lookup
 return gf2_64_basis::nim_to_poly(u64(inv_n16));
}
u64 inv(u64 a) {
 const u64 b1= frob16(a);
 const u64 b2= frob16(b1);
 const u64 b3= frob16(b2);
 const u64 beta= mul(mul(b1, b2), b3);
 const u64 N= mul(a, beta);
 return mul(beta, inv_in_f16_direct(N));
}
}  // namespace gf2_64_pclmul_norm_inv16
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_inv16::init_tables;
  using gf2_64_pclmul_norm_inv16::inv;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
