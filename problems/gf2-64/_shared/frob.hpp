#pragma once
// GF(2^64) の Frobenius^k (= a^{2^k}) を byte table で計算する current best 実装。
// 提供: frob4, frob8, frob16, frob32 (k = 2, 3, 4, 5)
//
// 各 frob_k は F_2-線型写像なので 8 byte tables (= 16 KB / k) で 8 lookup + 8 XOR。
// テーブルは constexpr で compile-time に作られる (PDEP は使わず bit-interleave で
// 構成して constexpr 評価可)。
//
// 利用側ルール:
//   gf2-64-frob{4,8,16,32}/algos/* は本ファイルを使ってはいけない (frob_k の
//     比較対象なので)
//   それ以外の problem (div / pow / sqrt / log) は building block として OK
//   k=2 (sq) は _shared/sq.hpp が PDEP で提供するのでここには無い
//
// 更新時の手順:
//   1. gf2-64-frob{k}/algos/ で新バリアントを追加し比較
//   2. 速い側に確定したら本ファイルを更新
#include "_common.hpp"
namespace gf2_64_pclmul {
namespace _frob_detail {
constexpr u8 RED_TABLE[]= {0, 27, 45, 54, 90, 65, 119, 108};
constexpr u64 spread_constexpr(u32 a) {
 u64 x= a;
 x= (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x= (x | (x << 8)) & 0x00FF00FF00FF00FFull;
 x= (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
 x= (x | (x << 2)) & 0x3333333333333333ull;
 x= (x | (x << 1)) & 0x5555555555555555ull;
 return x;
}
constexpr u64 sq_constexpr(u64 a) {
 u64 h= spread_constexpr(u32(a >> 32));
 u64 d= h ^ (h << 1);
 return spread_constexpr(u32(a)) ^ RED_TABLE[h >> 60] ^ d ^ (d << 3);
}
template<int sq_count>
constexpr array<array<u64, 256>, 8> make_frob_table() {
 array<array<u64, 256>, 8> t{};
 for(int p= 0; p < 8; ++p) {
  for(int b= 0; b < 256; ++b) {
   u64 v= u64(b) << (8 * p);
   for(int i= 0; i < sq_count; ++i) v= sq_constexpr(v);
   t[p][b]= v;
  }
 }
 return t;
}
}  // namespace _frob_detail

inline constexpr auto FROB4_BYTE = _frob_detail::make_frob_table<2>();
inline constexpr auto FROB8_BYTE = _frob_detail::make_frob_table<3>();
inline constexpr auto FROB16_BYTE= _frob_detail::make_frob_table<4>();
inline constexpr auto FROB32_BYTE= _frob_detail::make_frob_table<5>();

[[gnu::always_inline]] inline u64 apply_frob_byte(const array<array<u64, 256>, 8>& t, u64 a) {
 return t[0][u8(a)] ^ t[1][u8(a >> 8)] ^ t[2][u8(a >> 16)] ^ t[3][u8(a >> 24)] ^ t[4][u8(a >> 32)] ^ t[5][u8(a >> 40)] ^ t[6][u8(a >> 48)] ^ t[7][u8(a >> 56)];
}
inline u64 frob4 (u64 a) { return apply_frob_byte(FROB4_BYTE,  a); }
inline u64 frob8 (u64 a) { return apply_frob_byte(FROB8_BYTE,  a); }
inline u64 frob16(u64 a) { return apply_frob_byte(FROB16_BYTE, a); }
inline u64 frob32(u64 a) { return apply_frob_byte(FROB32_BYTE, a); }
}
