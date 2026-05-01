#pragma once
// algos 共通: typedef とよく使うヘッダ。
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// この問題は GF(2) (mod 2) 上の N×N {0,1} 行列の逆行列。1 ≤ N ≤ 4096。
// 逆行列が存在すれば N 行の '0'/'1' 列を返す。存在しなければ空 vector を返し、
// base.cpp 側で "-1" と出力する。
