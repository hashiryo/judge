#pragma once
// algos 共通: typedef とよく使うヘッダ。
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// この問題は mod 2 の {0,1} 行列の積。1 ≤ N, M, K ≤ 4096。
// 入力: N×M と M×K の行列 (各セルは '0' または '1')
// 出力: N×K の積行列
