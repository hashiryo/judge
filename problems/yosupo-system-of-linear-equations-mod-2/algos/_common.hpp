#pragma once
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// この問題は GF(2) 上の連立 1 次方程式 Ax = b。1 ≤ N, M ≤ 4096。
// 出力:
//   解なし: 空 vector<string> を返す → base.cpp で "-1" を出力
//   解あり: R+1 行 (R = null space の次元 = 自由変数の数)
//     [0] 行目: 特解
//     [1..R] 行目: null space の基底 (各行は M bit の '0'/'1' 文字列)
