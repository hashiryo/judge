#pragma once
// pclmul_norm_inv16.hpp の改善版:
// nim_to_poly は入力が「上位 48 bit = 0」と分かっているので 8 lookups のうち 2 個
// だけで済む。明示的に embed_f16 (2 lookups only) を書くことで確実に高速化。
// 同様に poly_to_nim は subfield 元なので低 16 bit しか使わない → extract_f16
// (8 lookups で 16-bit 出力) で 64-bit XOR を 16-bit XOR に圧縮。
//
// 期待: 6 lookups 節約 (per inv)、~18 cycle 改善。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
#include "../../_shared/basis_change.hpp"
namespace gf2_64_pclmul_norm_inv16_compact {
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
using gf2_64_pclmul::frob16;
// poly basis 64-bit (subfield 元) → 16-bit nim 表現
inline u16 EXTRACT_F16_BYTE[8][256];
// nim 16-bit → poly basis 64-bit
inline u64 EMBED_F16_BYTE[2][256];
// F_{2^16} 直接 inv テーブル (nim 基底)
inline u16 INV16[65536];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 // EXTRACT_F16 = low 16 bits of BASIS_BYTE
 for(int p= 0; p < 8; ++p) {
  for(int b= 0; b < 256; ++b) {
   EXTRACT_F16_BYTE[p][b]= u16(gf2_64_basis::BASIS_BYTE[p][b]);
  }
 }
 // EMBED_F16 = INV_BYTE[0..1] (低 2 byte に対する nim_to_poly 寄与)
 for(int p= 0; p < 2; ++p) {
  for(int b= 0; b < 256; ++b) {
   EMBED_F16_BYTE[p][b]= gf2_64_basis::INV_BYTE[p][b];
  }
 }
 // F_{2^16} の log/exp を一旦構築 → INV16 を作る
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
[[gnu::always_inline]] inline u16 extract_f16(u64 N) { return EXTRACT_F16_BYTE[0][u8(N)] ^ EXTRACT_F16_BYTE[1][u8(N >> 8)] ^ EXTRACT_F16_BYTE[2][u8(N >> 16)] ^ EXTRACT_F16_BYTE[3][u8(N >> 24)] ^ EXTRACT_F16_BYTE[4][u8(N >> 32)] ^ EXTRACT_F16_BYTE[5][u8(N >> 40)] ^ EXTRACT_F16_BYTE[6][u8(N >> 48)] ^ EXTRACT_F16_BYTE[7][u8(N >> 56)]; }
[[gnu::always_inline]] inline u64 embed_f16(u16 n) { return EMBED_F16_BYTE[0][u8(n)] ^ EMBED_F16_BYTE[1][u8(n >> 8)]; }
u64 inv_in_f16(u64 N_poly) {
 const u16 n16= extract_f16(N_poly);
 const u16 inv_n16= INV16[n16];
 return embed_f16(inv_n16);
}
u64 inv(u64 a) {
 const u64 b1= frob16(a);
 const u64 b2= frob16(b1);
 const u64 b3= frob16(b2);
 const u64 beta= mul(mul(b1, b2), b3);
 const u64 N= mul(a, beta);
 return mul(beta, inv_in_f16(N));
}
}  // namespace gf2_64_pclmul_norm_inv16_compact
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_inv16_compact::init_tables;
  using gf2_64_pclmul_norm_inv16_compact::inv;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
