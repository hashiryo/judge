#pragma once
// GF(2^64) の current best inv: Fermat 風 a^(2^64 - 2)。 pow 経由。
// (Itoh-Tsujii や norm-based 法のほうが速いはずだが、 比較は gf2-64-div/ で行うので
//  ここでは素朴 baseline のみ提供。 速い実装が確定したらここを差し替える。)
//
// 利用側ルール: gf2-64-div/algos/* は本ファイルを使ってはいけない (inv の比較対象なので)。
// それ以外の problem は building block として使用 OK。
#include "pow.hpp"
namespace gf2_64_pclmul {
[[gnu::target("pclmul")]] inline u64 inv(u64 a) { return pow(a, ~(u64)1); }
}
