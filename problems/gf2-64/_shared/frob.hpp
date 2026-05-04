#pragma once
// GF(2^64) の Frobenius^k (= a^{2^k}) を byte table で計算する current best 実装。
//
// 命名規約: frobK の K は **sq の適用回数** (= 結果は a^{2^K} に等しい)。
//   frob1 (= sq, _shared/sq.hpp で提供)
//   frob2 (= a^4)
//   frob4 (= a^16)
//   frob8 (= a^256)
//   frob16 (= a^65536)
//
// 各 frobK は F_2-線型写像なので 8 byte tables (= 16 KB / table) で 8 lookup + 7 XOR。
// テーブルは constexpr で compile-time に作られる (PDEP は使わず bit-interleave で
// 構成して constexpr 評価可)。
//
// 利用側ルール:
//   gf2-64-frob{2,4,8,16}/algos/* は本ファイルを使ってはいけない (frobK の
//     比較対象なので)。
//   それ以外の problem (div / pow / sqrt / log) は building block として OK。
//   K=1 (= sq) は _shared/sq.hpp が PDEP で提供。 ここには無い。
//
// Itoh-Tsujii で必要な K = 1, 2, 4, 8, 16 を全カバー。
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
template<int sq_count> constexpr array<array<u64, 256>, 8> make_frob_table() {
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

// FROBK_BYTE[p][b] = frobK(b << 8p)。 K = sq 適用回数。
inline constexpr auto FROB2_BYTE = _frob_detail::make_frob_table<2>();
inline constexpr auto FROB4_BYTE = _frob_detail::make_frob_table<4>();
inline constexpr auto FROB8_BYTE = _frob_detail::make_frob_table<8>();
inline constexpr auto FROB16_BYTE= _frob_detail::make_frob_table<16>();

// frob32, frob48 は norm decomposition で使う (α + α^{2^16} + α^{2^32} + α^{2^48})
// 直接 sq×32 / sq×48 で constexpr 評価すると step 上限超過する可能性があるので、
// 既存の FROB16_BYTE を chain で適用して構築 (各 entry 16-32 ops で軽い)。
namespace _frob_detail {
constexpr u64 apply_frob16_constexpr(u64 a) {
 return FROB16_BYTE[0][u8(a)] ^ FROB16_BYTE[1][u8(a >> 8)] ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)] ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)] ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}
constexpr array<array<u64, 256>, 8> make_frob_chain(int chain_count) {
 array<array<u64, 256>, 8> t{};
 for(int p= 0; p < 8; ++p) {
  for(int b= 0; b < 256; ++b) {
   u64 v= u64(b) << (8 * p);
   for(int i= 0; i < chain_count; ++i) v= apply_frob16_constexpr(v);
   t[p][b]= v;
  }
 }
 return t;
}
}  // namespace _frob_detail
inline constexpr auto FROB32_BYTE= _frob_detail::make_frob_chain(2);  // sq×32 = frob16 ∘ frob16
inline constexpr auto FROB48_BYTE= _frob_detail::make_frob_chain(3);  // sq×48 = frob16^3

// 各 frobK は対応する byte table を直接展開 (関数引数経由の indirection を避ける)。
inline u64 frob2(u64 a) {
 return FROB2_BYTE[0][u8(a)] ^ FROB2_BYTE[1][u8(a >> 8)] ^ FROB2_BYTE[2][u8(a >> 16)] ^ FROB2_BYTE[3][u8(a >> 24)] ^ FROB2_BYTE[4][u8(a >> 32)] ^ FROB2_BYTE[5][u8(a >> 40)] ^ FROB2_BYTE[6][u8(a >> 48)] ^ FROB2_BYTE[7][u8(a >> 56)];
}
inline u64 frob4(u64 a) {
 return FROB4_BYTE[0][u8(a)] ^ FROB4_BYTE[1][u8(a >> 8)] ^ FROB4_BYTE[2][u8(a >> 16)] ^ FROB4_BYTE[3][u8(a >> 24)] ^ FROB4_BYTE[4][u8(a >> 32)] ^ FROB4_BYTE[5][u8(a >> 40)] ^ FROB4_BYTE[6][u8(a >> 48)] ^ FROB4_BYTE[7][u8(a >> 56)];
}
inline u64 frob8(u64 a) {
 return FROB8_BYTE[0][u8(a)] ^ FROB8_BYTE[1][u8(a >> 8)] ^ FROB8_BYTE[2][u8(a >> 16)] ^ FROB8_BYTE[3][u8(a >> 24)] ^ FROB8_BYTE[4][u8(a >> 32)] ^ FROB8_BYTE[5][u8(a >> 40)] ^ FROB8_BYTE[6][u8(a >> 48)] ^ FROB8_BYTE[7][u8(a >> 56)];
}
inline u64 frob16(u64 a) {
 return FROB16_BYTE[0][u8(a)] ^ FROB16_BYTE[1][u8(a >> 8)] ^ FROB16_BYTE[2][u8(a >> 16)] ^ FROB16_BYTE[3][u8(a >> 24)] ^ FROB16_BYTE[4][u8(a >> 32)] ^ FROB16_BYTE[5][u8(a >> 40)] ^ FROB16_BYTE[6][u8(a >> 48)] ^ FROB16_BYTE[7][u8(a >> 56)];
}
inline u64 frob32(u64 a) {
 return FROB32_BYTE[0][u8(a)] ^ FROB32_BYTE[1][u8(a >> 8)] ^ FROB32_BYTE[2][u8(a >> 16)] ^ FROB32_BYTE[3][u8(a >> 24)] ^ FROB32_BYTE[4][u8(a >> 32)] ^ FROB32_BYTE[5][u8(a >> 40)] ^ FROB32_BYTE[6][u8(a >> 48)] ^ FROB32_BYTE[7][u8(a >> 56)];
}
inline u64 frob48(u64 a) {
 return FROB48_BYTE[0][u8(a)] ^ FROB48_BYTE[1][u8(a >> 8)] ^ FROB48_BYTE[2][u8(a >> 16)] ^ FROB48_BYTE[3][u8(a >> 24)] ^ FROB48_BYTE[4][u8(a >> 32)] ^ FROB48_BYTE[5][u8(a >> 40)] ^ FROB48_BYTE[6][u8(a >> 48)] ^ FROB48_BYTE[7][u8(a >> 56)];
}
}
