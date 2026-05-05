#pragma once
// Norm-based inv で F_{2^16} subfield を log/exp テーブル直引きする版。
//
// pclmul_norm.hpp は subfield inv に Itoh-Tsujii (6 muls + 15 sqs ≈ 60 cycles) を
// 使うが、F_{2^16} の log/exp テーブル (256 KiB) を持てば 2 lookup (~20 cycles) で済む。
//
// 鍵: nim 基底では F_{2^16} subfield = low 16 bits。既存の basis_change.hpp で
//     poly ↔ nim を 8 byte lookup × 2 で行えるので、subfield inv 全体は
//       poly_to_nim → low 16 bits → log/exp inv → upper 0 → nim_to_poly
//     合計 8 + 2 + 8 = 18 byte lookups (~50 cycles)。
//
// log/exp テーブルは Nimber.hpp と同じ recurrence で構築 (PW[i] = nim 表現での s^i)。

#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
#include "../../_shared/basis_change.hpp"
namespace gf2_64_pclmul_norm_logexp {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
// F_{2^16} log/exp テーブル (Nimber 互換, nim 基底での値が index/値となる)
inline u16 PW16[65536], LN16[65536];
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 // F_{2^16} log/exp テーブル (Nimber.hpp と同じ recurrence)。
 // 第 1 段で poly 基底 (s = 2 を生成元、r(s) = s^16 + ...) で PW[i] = s^i
 PW16[0]= PW16[65535]= 1;
 for(int i= 1; i < 65535; ++i) {
  PW16[i]= u16((PW16[i - 1] << 1) ^ (0x1681fu & u16(-(PW16[i - 1] >= 0x8000u))));
 }
 // 第 2 段で nim 基底に変換し、変換後の値で LN を引けるようにする
 constexpr u16 f2n[16]= {0x0001u, 0x2827u, 0x392bu, 0x8000u, 0x20fdu, 0x4d1du, 0xde4au, 0x0a17u, 0x3464u, 0xe3a9u, 0x6d8du, 0x34bcu, 0xa921u, 0xa173u, 0x0ebcu, 0x0e69u};
 for(int i= 1; i < 65535; ++i) {
  u16 x= PW16[i], y= 0;
  for(; x; x&= x - 1) y^= f2n[__builtin_ctz(x)];
  PW16[i]= y;
  LN16[y]= u16(i);
 }
 LN16[1]= 0;
}
// F_{2^16} subfield 元 N (poly 基底 64-bit) の inv を log/exp で計算 → 64-bit poly に戻す
[[gnu::target("pclmul")]] u64 inv_in_f16_logexp(u64 N_poly) {
 const u64 N_nim= gf2_64_basis::poly_to_nim(N_poly);
 const u16 n16= u16(N_nim);  // nim 基底では subfield 元の上位 48 bit は 0
 const u16 log_n= LN16[n16];
 const u16 inv_log= u16((65535u - u32(log_n)) % 65535u);
 const u16 inv_n16= PW16[inv_log];
 return gf2_64_basis::nim_to_poly(u64(inv_n16));
}
[[gnu::target("pclmul")]] u64 inv(u64 a) {
 const u64 b1= frob16(a);
 const u64 b2= frob16(b1);
 const u64 b3= frob16(b2);
 const u64 beta= mul(mul(b1, b2), b3);
 const u64 N= mul(a, beta);
 return mul(beta, inv_in_f16_logexp(N));
}
}  // namespace gf2_64_pclmul_norm_logexp
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  using gf2_64_pclmul_norm_logexp::init_tables;
  using gf2_64_pclmul_norm_logexp::inv;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
