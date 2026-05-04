#pragma once
// GF(2^64) の current best sq: PDEP (BMI2) で bit を spread → P(x) で reduce。
// GF(2) では (Σ a_i x^i)^2 = Σ a_i x^{2i} (cross term ゼロ) なので、 PDEP 1 命令で書ける。
// BMI2 が無い環境では bit interleave fallback。
//
// 利用側ルール: gf2-64-sq/algos/* / gf2-64-mul/algos/* は本ファイルを使ってはいけない。
// それ以外の problem (div/pow/sqrt/log) は building block として使用 OK。
#include "_common.hpp"
#include "mul.hpp"  // RED 定数を共有
namespace gf2_64_pclmul {
inline u64 spread_bits(u32 a) {
#ifdef __BMI2__
 return _pdep_u64(u64(a), 0x5555555555555555ull);
#else
 // fallback: bit interleave for 32-bit input
 u64 x= a;
 x= (x | (x << 16)) & 0x0000FFFF0000FFFFull;
 x= (x | (x << 8)) & 0x00FF00FF00FF00FFull;
 x= (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
 x= (x | (x << 2)) & 0x3333333333333333ull;
 x= (x | (x << 1)) & 0x5555555555555555ull;
 return x;
#endif
}
inline u64 sq(u64 a) {
 u64 h= spread_bits(a >> 32), d= h ^ (h << 1);
 return spread_bits(a) ^ RED[h >> 60] ^ d ^ (d << 3);
}
}
