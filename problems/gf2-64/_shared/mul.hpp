#pragma once
// GF(2^64) の current best mul: PCLMUL + 多項式 P(x) = x^64+x^4+x^3+x+1 で reduce。
//
// 利用側ルール: gf2-64-mul/algos/* は本ファイルを使ってはいけない (mul の比較対象なので)。
// それ以外の problem (sq/div/pow/sqrt/log) は building block として使用 OK。
#include "_common.hpp"
namespace gf2_64_pclmul {
constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
PCLMUL inline u64 mul(u64 a, u64 b) {
 __m128i v= _mm_clmulepi64_si128(_mm_cvtsi64_si128(a), _mm_cvtsi64_si128(b), 0);
 u64 h= (u64)v[1], d= h ^ (h << 1);
 return (u64)v[0] ^ RED[h >> 60] ^ d ^ (d << 3);
}
}
