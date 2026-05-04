#pragma once
// GF(2^64) の current best sqrt: a^(2^63) を pow 経由 (= 63 sq)。
// (Frobenius byte-table linear map のほうが速いはずだが、 比較は gf2-64-sqrt/ で行う。)
//
// 利用側ルール: gf2-64-sqrt/algos/* は本ファイルを使ってはいけない (sqrt の比較対象なので)。
// それ以外の problem は building block として使用 OK。
#include "pow.hpp"
namespace gf2_64_pclmul {
[[gnu::target("pclmul")]] inline u64 sqrt(u64 a) { return pow(a, (u64)1 << 63); }
}
