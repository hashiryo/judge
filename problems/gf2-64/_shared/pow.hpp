#pragma once
// GF(2^64) の current best pow: 二進展開 (mul + sq の current best を使う)。
//
// 利用側ルール: gf2-64-pow/algos/* は本ファイルを使ってはいけない (pow の比較対象なので)。
// それ以外の problem (log 等) は building block として使用 OK。
#include "mul.hpp"
#include "sq.hpp"
namespace gf2_64_pclmul {
[[gnu::target("pclmul")]] inline u64 pow(u64 a, u64 e) {
 u64 res= 1;
 while(e) {
  if(e & 1) res= mul(res, a);
  a= sq(a);
  e>>= 1;
 }
 return res;
}
}
